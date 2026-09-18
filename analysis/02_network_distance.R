# 02_network_distance.R
# Street-network distance from each tract's Census 2010 population-weighted center to the
# nearest recreational spaces (mean of the three nearest, and the single nearest), measured
# along the TIGER/Line 2016 road network of Los Angeles County, and the equivalent walking time.
#
# Inputs : ../public-sources/downloads/tl_2016_06037_edges.zip   (TIGER/Line 2016 All Lines edges; see MANIFEST.csv)
#          ../public-sources/downloads/parks_recreational_spaces.geojson  (the study park layer)
#          ../public-sources/downloads/pna2016_parks_open_space.geojson   (the county's 2016 park inventory; sensitivity analysis)
#          ../public-sources/derived/tract_exposures.csv          (population-weighted centers; built by 01_prepare_inputs.R)
# Output : ../public-sources/derived/tract_network_distance.csv
#          (the primary measure uses all 2,509 features of the study layer; sensitivity variants use open-access
#           features only, all features with overlapping polygons merged, and the 2016 county park inventory)
# Needs  : R with sf, and Python 3 (standard library only) for the shortest-path step (analysis/network_distance.py).
#
# Method. Road edges (ROADFLG = Y) other than primary roads (limited-access highways, MTFCC S1100) and
# ramps (S1630) form an undirected network whose nodes are the TIGER node identifiers and whose edge
# weights are edge lengths in California Albers (EPSG:3310). Only the connected parts of the network with
# at least 500 nodes are used (the main network, which holds 99.5% of the nodes, and Santa Catalina Island);
# smaller disconnected fragments are ignored. Park access points are every used node inside a park polygon
# or within 60 m of its boundary, every point where a used road edge crosses the park boundary, and the
# nearest point on the nearest used road edge; each tract's population-weighted center is joined to the
# nearest point on the nearest used road edge, and a center that itself lies inside a park or within 60 m
# of its boundary reaches that park directly. The straight-line gaps between a center and the road, and
# between the road and a park boundary, are added to the shortest path. Walking time assumes 3 miles per hour.
#
# Run from the repository root, after 01_prepare_inputs.R:
#   Rscript analysis/02_network_distance.R > analysis/output/02_network_distance.log 2>&1
# The intermediate network files (nodes, edges, access points) are written to a temporary directory;
# to keep them, give a directory as the first argument.
suppressPackageStartupMessages({ library(sf) })
options(width = 120, stringsAsFactors = FALSE)

here <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) dirname(normalizePath(f)) else getwd()
})
PKG <- normalizePath(file.path(here, ".."))
SRC <- file.path(PKG, "public-sources", "downloads")
DER <- file.path(PKG, "public-sources", "derived")
PY  <- file.path(here, "network_distance.py")

MI <- 1609.344                       # meters per mile
BAND_M <- 60                         # nodes within this distance of a park polygon are access points
MIN_NODES <- 500                     # connected parts of the network smaller than this are ignored
check <- function(ok, msg) {
  if (!isTRUE(ok)) stop("Check failed: ", msg, call. = FALSE)
  cat("  ok:", msg, "\n")
}
hr <- function(t) cat("\n== ", t, " ", strrep("=", max(3, 90 - nchar(t))), "\n", sep = "")
args <- commandArgs(TRUE)
tmp <- if (length(args)) args[1] else file.path(tempdir(), "network")
dir.create(tmp, showWarnings = FALSE, recursive = TRUE)

hr("1. Road network (TIGER/Line 2016 All Lines edges, Los Angeles County)")
ed <- file.path(tempdir(), "edges"); dir.create(ed, showWarnings = FALSE)
unzip(file.path(SRC, "tl_2016_06037_edges.zip"), exdir = ed)
e <- st_read(list.files(ed, pattern = "[.]shp$", full.names = TRUE)[1], quiet = TRUE)
cat("  edges in the file:", nrow(e), "| road edges (ROADFLG = Y):", sum(e$ROADFLG == "Y"), "\n")
rd <- e[e$ROADFLG == "Y", ]
cat("  road edges by MTFCC feature class:\n"); print(table(rd$MTFCC))
DROP <- c("S1100", "S1630")           # primary roads (limited-access highways) and ramps
w <- rd[!(rd$MTFCC %in% DROP), ]
cat("  dropped:", sum(rd$MTFCC == "S1100"), "primary-road edges and", sum(rd$MTFCC == "S1630"), "ramp edges | kept:", nrow(w), "\n")
w <- st_transform(w, 3310)
check(all(st_geometry_type(w) == "LINESTRING"), "every road edge is a single linestring")
w$len_m <- as.numeric(st_length(w))
cat("  network length:", round(sum(w$len_m) / MI), "miles | zero-length edges:", sum(w$len_m == 0), "\n")

# Nodes: TIGER node identifiers at the first (TNIDF) and last (TNIDT) vertex of each edge
co <- st_coordinates(st_geometry(w))
first <- !duplicated(co[, "L1"]); last <- !duplicated(co[, "L1"], fromLast = TRUE)
tnidf <- as.numeric(w$TNIDF); tnidt <- as.numeric(w$TNIDT)
ends <- rbind(data.frame(tnid = tnidf, x = co[first, "X"], y = co[first, "Y"]),
              data.frame(tnid = tnidt, x = co[last, "X"], y = co[last, "Y"]))
spread <- tapply(seq_len(nrow(ends)), ends$tnid, function(i) if (length(i) > 1) max(dist(ends[i, c("x", "y")])) else 0)
check(max(spread) < 0.01, "each TIGER node identifier has a single coordinate")
nd <- ends[!duplicated(ends$tnid), ]
nd$node_id <- seq_len(nrow(nd)) - 1L
w$u <- nd$node_id[match(tnidf, nd$tnid)]; w$v <- nd$node_id[match(tnidt, nd$tnid)]
cat("  nodes:", nrow(nd), "| edges:", nrow(w), "\n")

# Connected parts of the network (label propagation with pointer jumping)
comp_of <- local({
  p <- seq_len(nrow(nd)); u <- w$u + 1L; v <- w$v + 1L
  repeat {
    ru <- p[u]; rv <- p[v]; d <- ru != rv
    if (!any(d)) break
    hi <- pmax(ru[d], rv[d]); lo <- pmin(ru[d], rv[d])
    m <- tapply(lo, hi, min); idx <- as.integer(names(m))
    p[idx] <- pmin(p[idx], as.integer(m))
    repeat { pp <- p[p]; if (identical(pp, p)) break; p <- pp }
  }
  p
})
csize <- table(comp_of)
big <- as.integer(names(csize)[csize >= MIN_NODES])
cat("  connected parts:", length(csize), "| largest:", max(csize), "nodes (", round(100 * max(csize) / nrow(nd), 2),
    "% ) | parts with at least", MIN_NODES, "nodes:", length(big), "( sizes", paste(sort(csize[csize >= MIN_NODES], decreasing = TRUE), collapse = ", "),
    ") | next largest:", paste(head(sort(csize[csize < MIN_NODES], decreasing = TRUE), 5), collapse = ", "), "\n")
nd$component <- match(comp_of, big)                  # 1 = main network, 2 = the next used part, NA = ignored fragment
w$component <- nd$component[w$u + 1L]
used_n <- !is.na(nd$component); used_e <- !is.na(w$component)
cat("  used:", sum(used_n), "nodes and", sum(used_e), "edges; ignored fragments:", sum(!used_n), "nodes,", sum(!used_e), "edges\n")
nds <- st_as_sf(nd[used_n, ], coords = c("x", "y"), crs = 3310)
wr <- w[used_e, ]
write.csv(data.frame(node_id = nds$node_id, component = nds$component), file.path(tmp, "nodes.csv"), row.names = FALSE)
write.csv(data.frame(edge_id = seq_len(nrow(wr)) - 1L, u = wr$u, v = wr$v, len_m = round(wr$len_m, 3)), file.path(tmp, "edges.csv"), row.names = FALSE)

# nearest used edge for a set of geometries: straight-line gap and position of the nearest point along the edge
snap_to_edge <- function(geom) {
  j <- st_nearest_feature(geom, st_geometry(wr))
  gap <- as.numeric(st_distance(geom, st_geometry(wr)[j], by_element = TRUE))
  np <- st_nearest_points(geom, st_geometry(wr)[j], pairwise = TRUE)
  on_edge <- st_cast(np, "POINT")[2 * seq_along(j)]              # the end of each connecting line that lies on the edge
  pos <- as.numeric(st_line_project(st_geometry(wr)[j], on_edge))
  data.frame(edge_id = j - 1L, gap_m = gap, pos_m = pmin(pmax(pos, 0), wr$len_m[j]), component = wr$component[j])
}

hr("2. Park access points (the study park layer, 2,509 polygons)")
pk <- st_read(file.path(SRC, "parks_recreational_spaces.geojson"), quiet = TRUE)
pk <- st_make_valid(st_transform(pk, 3310))
check(nrow(pk) == 2509, "2,509 recreational-space polygons")
cat("  access type of the features:\n"); print(table(pk$ACCESS_TYP, useNA = "ifany"))
# Park sets from the study layer: all features (the primary measure, as in the straight-line measure); open-access features only;
# all features with overlapping polygons merged, so that duplicate polygons of one park from different
# source datasets cannot count as two of the three nearest parks
# merged set: features whose interiors overlap another feature (duplicate or partly duplicate polygons from
# different source datasets; features that merely touch are not merged) are grouped by the connected components
# of the overlap graph, and each group is unioned into one park; features that stand alone are kept as they are
grp <- local({
  nb <- st_intersects(st_geometry(pk)); tc <- st_touches(st_geometry(pk))
  nb <- lapply(seq_along(nb), function(i) setdiff(nb[[i]], c(i, tc[[i]])))
  p <- seq_len(nrow(pk))
  find <- function(i) { while (p[i] != i) { p[i] <<- p[p[i]]; i <- p[i] }; i }
  for (i in seq_along(nb)) for (j in nb[[i]]) if (j > i) { a <- find(i); b <- find(j); if (a != b) p[max(a, b)] <- min(a, b) }
  vapply(seq_along(nb), find, integer(1))
})
merged <- st_sfc(lapply(split(seq_len(nrow(pk)), grp), function(i) if (length(i) == 1) st_geometry(pk)[[i]] else st_union(st_geometry(pk)[i])[[1]]), crs = 3310)
# Fourth park set: the county's 2016 Parks and Recreation Needs Assessment inventory (all features except
# 'No Public Access', as in 01_prepare_inputs.R), so that parks, roads (2016) and outcomes (2015-2017) share one period
p16 <- st_read(file.path(SRC, "pna2016_parks_open_space.geojson"), quiet = TRUE)
p16 <- st_make_valid(st_transform(p16[!(p16$ACCESS_TYP %in% "No Public Access"), ], 3310))
sets <- list(all = st_geometry(pk),
             open = st_geometry(pk)[pk$ACCESS_TYP %in% "Open Access"],
             merged = st_make_valid(merged),
             pna2016 = st_geometry(p16))
gsize <- table(grp)
cat("  park polygons per set: all", length(sets$all), "| open access", length(sets$open), "| merged (overlapping features grouped)",
    length(sets$merged), "| features in a group of two or more:", sum(gsize[as.character(grp)] > 1), "| largest group:", max(gsize), "\n")
check(length(sets$open) == 2230 && length(sets$merged) < 2509, "open-access subset has 2,230 features; merging reduces the count")
cat("  2016 county park inventory: features used (all except 'No Public Access'):", length(sets$pna2016), "\n")
access_points <- function(g, tag) {
  hits <- st_intersects(st_buffer(g, BAND_M), nds)
  pn <- data.frame(park_id = rep(seq_along(g), lengths(hits)), node_id = nds$node_id[unlist(hits)])
  pn$gap_m <- unlist(lapply(seq_along(g), function(i) if (length(hits[[i]])) as.numeric(st_distance(g[i], st_geometry(nds)[hits[[i]]])) else numeric(0)))
  pe <- snap_to_edge(g)
  # every point where a used road edge crosses the park boundary is an access point with no gap
  bd <- st_boundary(g); cx <- st_intersects(bd, st_geometry(wr))
  pairs <- data.frame(park_id = rep(seq_along(g), lengths(cx)), edge_id = unlist(cx))
  xy <- lapply(seq_len(nrow(pairs)), function(k) {          # crossing points of each boundary-edge pair
    x <- st_intersection(bd[[pairs$park_id[k]]], st_geometry(wr)[[pairs$edge_id[k]]])
    if (st_is_empty(x)) matrix(numeric(0), 0, 2) else unname(st_coordinates(x)[, 1:2, drop = FALSE])
  })
  npt <- vapply(xy, nrow, integer(1))
  pc_pts <- st_sfc(lapply(seq_len(sum(npt)), function(j) NULL), crs = 3310)
  pc_pts <- st_sfc(apply(do.call(rbind, xy), 1, st_point, simplify = FALSE), crs = 3310)
  cr <- data.frame(park_id = rep(pairs$park_id, npt), edge_id = rep(pairs$edge_id, npt))
  cr$pos_m <- as.numeric(st_line_project(st_geometry(wr)[cr$edge_id], pc_pts))
  cr$pos_m <- pmin(pmax(cr$pos_m, 0), wr$len_m[cr$edge_id])
  cat("  [", tag, "] used nodes inside or within ", BAND_M, " m of a park: ", nrow(pn), " pairs | parks with none: ", sum(lengths(hits) == 0),
      " | road crossings of park boundaries: ", nrow(cr), " points on ", length(unique(cr$park_id)), " parks",
      " | nearest used road edge per park: gap (m) median ", round(median(pe$gap_m), 1), ", 95th percentile ", round(quantile(pe$gap_m, .95), 1),
      ", max ", round(max(pe$gap_m), 1), " | parks by network part: ", paste(names(table(pe$component)), table(pe$component), sep = "=", collapse = ", "), "\n", sep = "")
  check(all(table(pe$component) >= 3), paste0("[", tag, "] every used part of the network reaches at least three parks"))
  write.csv(data.frame(park_id = pn$park_id, node_id = pn$node_id, gap_m = round(pn$gap_m, 3)), file.path(tmp, paste0("park_nodes_", tag, ".csv")), row.names = FALSE)
  write.csv(data.frame(park_id = c(seq_along(g), cr$park_id), edge_id = c(pe$edge_id, cr$edge_id - 1L), gap_m = round(c(pe$gap_m, rep(0, nrow(cr))), 3),
                       pos_m = round(c(pe$pos_m, cr$pos_m), 3)),
            file.path(tmp, paste0("park_edges_", tag, ".csv")), row.names = FALSE)
  # a tract center inside a park or within the band reaches that park directly (the same allowance as a node)
  op <- st_intersects(st_buffer(g, BAND_M), st_geometry(pc))
  od <- data.frame(park_id = rep(seq_along(g), lengths(op)), origin = unlist(op))
  od$gap_m <- unlist(lapply(seq_along(g), function(i) if (length(op[[i]])) as.numeric(st_distance(g[i], st_geometry(pc)[op[[i]]])) else numeric(0)))
  cat("  [", tag, "] tract centers inside or within ", BAND_M, " m of a park: ", length(unique(od$origin)), " (", sum(od$gap_m == 0), " center-park pairs inside)\n", sep = "")
  write.csv(data.frame(GEOID = ex$GEOID[od$origin], park_id = od$park_id, gap_m = round(od$gap_m, 3)), file.path(tmp, paste0("origin_parks_", tag, ".csv")), row.names = FALSE)
}
ex <- read.csv(file.path(DER, "tract_exposures.csv"), colClasses = c(GEOID = "character"))
check(nrow(ex) == 2343, "2,343 tracts in tract_exposures.csv")
pc <- st_as_sf(ex[, c("GEOID", "popcentroid_x", "popcentroid_y")], coords = c("popcentroid_x", "popcentroid_y"), crs = 3310)
for (tag in names(sets)) access_points(sets[[tag]], tag)

hr("3. Tract origins: Census 2010 population-weighted centers")
oe <- snap_to_edge(st_geometry(pc))
cat("  nearest used road edge per tract center: gap (m) median", round(median(oe$gap_m), 1), "| 95th percentile",
    round(quantile(oe$gap_m, .95), 1), "| max", round(max(oe$gap_m), 1), "(tract", ex$GEOID[which.max(oe$gap_m)], ")\n")
cat("  tract centers by network part:\n"); print(table(oe$component))
write.csv(data.frame(GEOID = ex$GEOID, edge_id = oe$edge_id, gap_m = round(oe$gap_m, 3), pos_m = round(oe$pos_m, 3)),
          file.path(tmp, "origin_edges.csv"), row.names = FALSE)

hr("4. Shortest paths along the network (Python 3, Dijkstra's algorithm)")
check(file.exists(PY), "analysis/network_distance.py is present")
res <- list()
for (tag in names(sets)) {
  status <- system2("python3", c(shQuote(PY), shQuote(tmp), tag), stdout = "", stderr = "")
  check(identical(as.integer(status), 0L), paste0("[", tag, "] network_distance.py finished without error"))
  r <- read.csv(file.path(tmp, paste0("network_result_", tag, ".csv")), colClasses = c(GEOID = "character", parks = "character"))
  check(nrow(r) == 2343 && identical(r$GEOID, ex$GEOID), paste0("[", tag, "] one result row per tract, in the same order"))
  res[[tag]] <- r
}

hr("5. Checks and output")
out <- data.frame(GEOID = res$all$GEOID,
                  net3_popcentroid = res$all$net3_m / MI, net1_popcentroid = res$all$net1_m / MI,
                  walk3_popcentroid_min = res$all$net3_m / MI * 60 / 3, walk1_popcentroid_min = res$all$net1_m / MI * 60 / 3,
                  net3_popcentroid_open = res$open$net3_m / MI, net3_popcentroid_merged = res$merged$net3_m / MI,
                  net3_popcentroid_pna2016 = res$pna2016$net3_m / MI,
                  net_snap_m = oe$gap_m, net_component = oe$component, net3_parks = res$all$parks)
check(!anyNA(out$net3_popcentroid) && !anyNA(out$net1_popcentroid) && !anyNA(out$net3_popcentroid_open) && !anyNA(out$net3_popcentroid_merged) &&
      !anyNA(out$net3_popcentroid_pna2016),
      "a network distance for every tract in every park set")
check(all(out$net3_popcentroid_open >= out$net3_popcentroid - 1e-6) && all(out$net3_popcentroid_merged >= out$net3_popcentroid - 1e-6),
      "restricting or merging the park set never shortens the three-park distance")
cat(sprintf("  three-park network distance, open-access parks only: median %.3f miles; r with the all-features measure %.3f\n",
            median(out$net3_popcentroid_open), cor(out$net3_popcentroid_open, out$net3_popcentroid)))
cat(sprintf("  three-park network distance, overlapping polygons merged: median %.3f miles; r with the all-features measure %.3f; tracts whose value changes by more than 0.01 mile: %d\n",
            median(out$net3_popcentroid_merged), cor(out$net3_popcentroid_merged, out$net3_popcentroid), sum(abs(out$net3_popcentroid_merged - out$net3_popcentroid) > 0.01)))
cat(sprintf("  three-park network distance to the 2016 county park inventory: median %.3f miles; r with the study-layer network measure %.3f; r with the straight-line distance to the same inventory from the tract centroid %.3f\n",
            median(out$net3_popcentroid_pna2016), cor(out$net3_popcentroid_pna2016, out$net3_popcentroid), cor(out$net3_popcentroid_pna2016, ex$dist3_pna2016)))
check(all(out$net1_popcentroid >= ex$dist1_popcentroid - 1e-6) && all(out$net3_popcentroid >= ex$dist3_popcentroid - 1e-6),
      "network distances are never shorter than the straight-line distances from the same centers")
ratio <- out$net3_popcentroid / ex$dist3_popcentroid
cat("  network / straight-line ratio (three-park mean, population-weighted center): quantiles\n")
print(round(quantile(ratio, c(0, .05, .25, .5, .75, .95, 1)), 3))
cat("  three-park network distance (miles): quantiles\n")
print(round(quantile(out$net3_popcentroid, c(0, .05, .25, .5, .75, .95, 1)), 3))
cat("  walking time to the three nearest parks at 3 mph (minutes): median", round(median(out$walk3_popcentroid_min), 1),
    "| quartiles", round(quantile(out$walk3_popcentroid_min, .25), 1), "-", round(quantile(out$walk3_popcentroid_min, .75), 1), "\n")
cat(sprintf("  Pearson r with the study exposure (geographic centroid, straight line): %.3f; with the straight-line distance from the same centers: %.3f\n",
            cor(out$net3_popcentroid, ex$recorded_dist3_miles), cor(out$net3_popcentroid, ex$dist3_popcentroid)))
cat("  tract centers joined to the road network by more than 0.1 mile:", sum(out$net_snap_m > 0.1 * MI),
    "| more than 0.25 mile:", sum(out$net_snap_m > 0.25 * MI), "\n")
write.csv(out, file.path(DER, "tract_network_distance.csv"), row.names = FALSE)
cat("  written:", file.path("public-sources", "derived", "tract_network_distance.csv"), "\n")

hr("Done")
print(sessionInfo())
