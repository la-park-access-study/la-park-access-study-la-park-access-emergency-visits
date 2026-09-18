# 01_prepare_inputs.R
# Builds the tract- and ZIP-level input tables used by analysis.R from
# the files in ../data (the author's original analysis files) and
# ../public-sources/downloads (public source data; see ../public-sources/MANIFEST.csv).
#
# Outputs (../public-sources/derived/):
#   ces4_official_la.csv       official CalEnviroScreen 4.0 values for the
#                              2,343 LA County tracts, including the age
#                              percentages (missing = NA)
#   ces3_ozone_la.csv          CalEnviroScreen 3.0 ozone (2012-2014) by tract
#   tract_exposures.csv        park-access metrics recomputed from geometry
#   zcta_tract_rel_la.csv      2010 ZCTA-to-tract population overlaps (LA)
#   cdph_asthma_zip_la.csv     CDPH asthma ED visit counts, LA ZIPs, 2015-2017
#   zcta_la_2010.gpkg          2010 ZCTA polygons for the LA ZCTAs
#   acs_uninsured_2015_2019.csv percentage uninsured by tract (ACS 2015-2019)
#
# Run from the repository root, with the MANIFEST files in public-sources/downloads/:
#   Rscript analysis/01_prepare_inputs.R > analysis/output/01_prepare_inputs.log 2>&1
suppressPackageStartupMessages({ library(sf); library(readxl) })
options(width = 120, stringsAsFactors = FALSE)

here <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) dirname(normalizePath(f)) else getwd()
})
PKG <- normalizePath(file.path(here, ".."))
DATA <- file.path(PKG, "data")
SRC <- file.path(PKG, "public-sources", "downloads")
DER <- file.path(PKG, "public-sources", "derived")
dir.create(DER, showWarnings = FALSE, recursive = TRUE)

MI <- 1609.344                       # meters per mile
check <- function(ok, msg) {
  if (!isTRUE(ok)) stop("Check failed: ", msg, call. = FALSE)
  cat("  ok:", msg, "\n")
}
geoid <- function(x) sprintf("%011.0f", round(as.numeric(x)))
mean3 <- function(r) mean(sort(r, partial = 1:3)[1:3])
hr <- function(t) cat("\n== ", t, " ", strrep("=", max(3, 90 - nchar(t))), "\n", sep = "")

hr("1. Tract layer (CalEnviroScreen 4.0 geometry, 2010 tracts)")
td <- file.path(tempdir(), "tracts"); dir.create(td, showWarnings = FALSE)
unzip(file.path(DATA, "revisedLAcensustracts.zip"), exdir = td)
tr <- st_read(list.files(td, pattern = "[.]shp$", full.names = TRUE)[1], quiet = TRUE)
tr <- st_make_valid(st_transform(tr[, "Tract"], 3310))
tr$GEOID <- geoid(tr$Tract)
check(nrow(tr) == 2343 && !anyDuplicated(tr$GEOID), "2,343 unique tracts")

# Centroids as ArcGIS Find Centroids ("contained within each input feature"):
# the geometric centroid, replaced by a point on the surface when it falls outside.
cen <- st_centroid(st_geometry(tr))
own <- mapply(function(w, i) i %in% w, st_within(cen, tr), seq_len(nrow(tr)))
cen[!own] <- st_point_on_surface(st_geometry(tr)[!own])
cat("  centroids replaced by point-on-surface:", sum(!own), "\n")

rec <- read.csv(file.path(DATA, "oct17csv.csv"), fileEncoding = "UTF-8-BOM")
rec$GEOID <- geoid(rec$Tract)
tr$recorded <- rec$MEAN_Total_Miles[match(tr$GEOID, rec$GEOID)]
check(!anyNA(tr$recorded), "recorded distances joined 1:1")

hr("2. Recreational Spaces layer (LA County DPH; the layer used in the study)")
pk <- st_read(file.path(SRC, "parks_recreational_spaces.geojson"), quiet = TRUE)
pk <- st_make_valid(st_transform(pk, 3310))
check(nrow(pk) == 2509, "2,509 recreational-space polygons")
print(table(pk$ACCESS_TYP, useNA = "ifany"))
a_pk <- as.numeric(st_area(pk))
cat("  zero-area polygons:", sum(a_pk == 0), "| exact duplicate geometries:", sum(duplicated(st_as_text(st_geometry(pk)))),
    "| features overlapping at least one other feature:", sum(lengths(st_overlaps(pk)) > 0), "\n")

D <- matrix(as.numeric(st_distance(cen, pk)), nrow = nrow(tr)) / MI
tr$dist3_planar <- apply(D, 1, mean3)
tr$dist1_planar <- apply(D, 1, min)

# Geodesic (spherical, s2) distances to the 8 nearest candidates per tract
sf_use_s2(TRUE)
cen_ll <- st_transform(cen, 4326)
pk_ll  <- st_make_valid(st_transform(st_geometry(pk), 4326))
cand <- t(apply(D, 1, function(r) order(r)[1:8]))
dg <- t(vapply(seq_len(nrow(tr)), function(i)
  as.numeric(st_distance(cen_ll[i], pk_ll[cand[i, ]])) / MI, numeric(8)))
tr$dist3_geodesic <- apply(dg, 1, mean3)

agree <- function(x, y) c(median_abs = median(abs(x - y)), mean_abs = mean(abs(x - y)),
                          r = cor(x, y), within_0.01 = mean(abs(x - y) <= 0.01),
                          within_0.02 = mean(abs(x - y) <= 0.02))
cmp <- rbind(planar_CA_Albers = agree(tr$dist3_planar, tr$recorded),
             geodesic_s2      = agree(tr$dist3_geodesic, tr$recorded))
cat("\nAgreement of recomputed 3-park mean distance with the recorded study values:\n")
print(round(cmp, 4))
check(max(cmp[, "r"]) >= 0.99, "recomputed distances reproduce the recorded exposure (r >= 0.99)")
# The same comparison with distances to park centroids instead of park edges
Dc <- matrix(as.numeric(st_distance(cen, st_centroid(st_geometry(pk)))), nrow = nrow(tr)) / MI
cat("\nAgreement if distances are measured to park centroids instead of park edges:\n")
print(round(agree(apply(Dc, 1, mean3), tr$recorded), 4))

hr("3. Population-weighted centroids (Census 2010 centers of population)")
cp <- read.csv(file.path(SRC, "cenpop2010_tr06.txt"),
               colClasses = c(STATEFP = "character", COUNTYFP = "character", TRACTCE = "character"))
cp <- cp[cp$COUNTYFP == "037", ]
cat("  2010 census tracts in Los Angeles County (Census centers-of-population file):", nrow(cp), "\n")
cp$GEOID <- paste0(cp$STATEFP, cp$COUNTYFP, cp$TRACTCE)
m <- match(tr$GEOID, cp$GEOID)
check(!anyNA(m), "every study tract has a 2010 center of population")
pc <- st_transform(st_as_sf(cp[m, ], coords = c("LONGITUDE", "LATITUDE"), crs = 4269), 3310)
Dp <- matrix(as.numeric(st_distance(st_geometry(pc), pk)), nrow = nrow(tr)) / MI
tr$dist3_popcentroid <- apply(Dp, 1, mean3)
tr$dist1_popcentroid <- apply(Dp, 1, min)
cat("  centroid-to-population-center shift (miles): median",
    round(median(as.numeric(st_distance(cen, st_geometry(pc), by_element = TRUE)) / MI), 3), "\n")

hr("4. Park area within 0.5 and 1 mile of the tract centroid")
pu <- st_union(st_geometry(pk))
for (r in c(0.5, 1)) {
  buf <- st_sf(id = seq_len(nrow(tr)), geometry = st_buffer(cen, r * MI))
  ix <- suppressWarnings(st_intersection(buf, st_sf(geometry = pu)))
  a <- tapply(as.numeric(st_area(ix)), ix$id, sum) / 4046.8564224   # m2 -> acres
  v <- numeric(nrow(tr)); v[as.integer(names(a))] <- a
  tr[[paste0("acres_", sub("[.]", "p", r), "mi")]] <- v
}

hr("5. 2016 Park Needs Assessment inventory and half-mile walksheds")
p16 <- st_read(file.path(SRC, "pna2016_parks_open_space.geojson"), quiet = TRUE)
print(table(p16$ACCESS_TYP, useNA = "ifany"))
p16 <- st_make_valid(st_transform(p16[!(p16$ACCESS_TYP %in% "No Public Access"), ], 3310))
cat("  2016 features used (all except 'No Public Access'):", nrow(p16), "\n")
D16 <- matrix(as.numeric(st_distance(cen, p16)), nrow = nrow(tr)) / MI
tr$dist3_pna2016 <- apply(D16, 1, mean3)

ws <- st_make_valid(st_transform(st_read(file.path(SRC, "pna2016_walkshed_halfmile.geojson"), quiet = TRUE), 3310))
wsu <- st_union(st_geometry(ws))
ix <- suppressWarnings(st_intersection(st_sf(id = seq_len(nrow(tr)), geometry = st_geometry(tr)),
                                       st_sf(geometry = wsu)))
a <- tapply(as.numeric(st_area(ix)), ix$id, sum)
share <- numeric(nrow(tr)); share[as.integer(names(a))] <- a
tr$walkshed2016_share <- pmin(1, share / as.numeric(st_area(tr)))
tr$popcentroid_in_walkshed2016 <- as.integer(lengths(st_intersects(st_geometry(pc), wsu)) > 0)

hr("6. Write tract exposure table")
xy  <- st_coordinates(cen); pxy <- st_coordinates(pc)
out <- data.frame(GEOID = tr$GEOID, recorded_dist3_miles = tr$recorded,
                  centroid_x = xy[, 1], centroid_y = xy[, 2], popcentroid_x = pxy[, 1], popcentroid_y = pxy[, 2],
                  dist3_planar = tr$dist3_planar, dist3_geodesic = tr$dist3_geodesic,
                  dist1_planar = tr$dist1_planar, dist3_popcentroid = tr$dist3_popcentroid,
                  dist1_popcentroid = tr$dist1_popcentroid,
                  acres_0p5mi = tr$acres_0p5mi, acres_1mi = tr$acres_1mi,
                  dist3_pna2016 = tr$dist3_pna2016, walkshed2016_share = tr$walkshed2016_share,
                  popcentroid_in_walkshed2016 = tr$popcentroid_in_walkshed2016)
write.csv(out, file.path(DER, "tract_exposures.csv"), row.names = FALSE)
print(summary(out[, -1]))

hr("7. Official CalEnviroScreen 4.0 and 3.0 values")
c4 <- suppressWarnings(read_excel(file.path(SRC, "ces40_results_datadictionary.xlsx"), sheet = "CES4.0FINAL_results"))
c4 <- c4[c4[["California County"]] %in% "Los Angeles", ]
c4o <- data.frame(GEOID = geoid(c4[["Census Tract"]]), population = c4[["Total Population"]],
                  asthma = c4[["Asthma"]], ami = c4[["Cardiovascular Disease"]],
                  poverty = c4[["Poverty"]], pm25 = c4[["PM2.5"]], ozone = c4[["Ozone"]])
check(nrow(c4o) == 2343 && all(c4o$GEOID %in% tr$GEOID), "2,343 official CES 4.0 LA tracts match the tract layer")
cat("  official missing values: asthma", sum(is.na(c4o$asthma)), "| AMI", sum(is.na(c4o$ami)),
    "| poverty", sum(is.na(c4o$poverty)), "| PM2.5", sum(is.na(c4o$pm25)), "| ozone", sum(is.na(c4o$ozone)), "\n")

# Age percentages come from the Demographic Profile sheet of the same workbook
num <- function(x) suppressWarnings(as.numeric(x))
dp <- suppressWarnings(read_excel(file.path(SRC, "ces40_results_datadictionary.xlsx"), sheet = "Demographic Profile", skip = 1))
dp <- dp[dp[["California County"]] %in% "Los Angeles", ]
dp <- data.frame(GEOID = geoid(dp[["Census Tract"]]), population_dp = num(dp[["Total Population"]]),
                 child_10 = num(dp[["Children < 10 years (%)"]]), elderly_65 = num(dp[["Elderly > 64 years (%)"]]))
c4o <- merge(c4o, dp, by = "GEOID", sort = FALSE)
check(nrow(c4o) == 2343 && all(c4o$population == c4o$population_dp), "Demographic Profile sheet: 2,343 tracts, same populations")
c4o$population_dp <- NULL
cat("  official missing age percentages:", sum(is.na(c4o$child_10)), "tracts, all without residents:",
    all(c4o$population[is.na(c4o$child_10)] == 0), "\n")
write.csv(c4o, file.path(DER, "ces4_official_la.csv"), row.names = FALSE)

c3 <- read.csv(file.path(SRC, "ces30_results_june2018.csv"), check.names = FALSE)
c3 <- c3[c3[["California County"]] == "Los Angeles", ]
c3o <- data.frame(GEOID = geoid(c3[["Census Tract"]]), ozone_ces3 = as.numeric(c3[["Ozone"]]))
check(nrow(c3o) == 2343 && all(tr$GEOID %in% c3o$GEOID), "CES 3.0 ozone matched 2,343/2,343")
write.csv(c3o, file.path(DER, "ces3_ozone_la.csv"), row.names = FALSE)

hr("8. ZIP-level inputs (CDPH asthma ED counts; 2010 ZCTA relationships and polygons)")
rel <- read.csv(file.path(SRC, "zcta_tract_rel_10.txt"), colClasses = c(ZCTA5 = "character", GEOID = "character"))
rel_la <- rel[rel$STATE == 6 & rel$COUNTY == 37, c("ZCTA5", "GEOID", "POPPT", "ZPOP", "ZPOPPCT")]
notin <- setdiff(tr$GEOID, rel_la$GEOID)
cat("  tracts with no ZCTA (no residents, so not in the primary sample):", paste(notin, collapse = ", "), "\n")
check(all(c4o$GEOID[c4o$population > 0] %in% rel_la$GEOID), "every tract with residents appears in the ZCTA relationship file")
write.csv(rel_la, file.path(DER, "zcta_tract_rel_la.csv"), row.names = FALSE)

az <- suppressWarnings(read_excel(file.path(SRC, "cdph_asthma_zip_allages.xls"), sheet = "All Ages with county"))
az <- az[az$County %in% "Los Angeles" & az$Year %in% 2015:2017, ]
zz <- aggregate(cbind(visits = Number_of_Asthma_ED_Visits, rate = Age_Adjusted_Rate_of_Asthma_ED_V) ~ Zip_Code,
                data = az, FUN = function(v) c(sum = sum(v), mean = mean(v), n = length(v)))
zip <- data.frame(zip = as.character(zz$Zip_Code), years = zz$visits[, "n"],
                  visits_2015_2017 = zz$visits[, "sum"], mean_age_adj_rate = zz$rate[, "mean"])
cat("  LA ZIPs 2015-2017:", nrow(zip), "| with all three years:", sum(zip$years == 3), "\n")
write.csv(zip, file.path(DER, "cdph_asthma_zip_la.csv"), row.names = FALSE)

zd <- file.path(tempdir(), "zcta"); dir.create(zd, showWarnings = FALSE)
unzip(file.path(SRC, "tl_2010_06_zcta510.zip"), exdir = zd)
zc <- st_read(list.files(zd, pattern = "[.]shp$", full.names = TRUE)[1], quiet = TRUE)
zc <- st_transform(zc[zc$ZCTA5CE10 %in% unique(rel_la$ZCTA5), "ZCTA5CE10"], 3310)
names(zc)[1] <- "ZCTA5"
if (file.exists(file.path(DER, "zcta_la_2010.gpkg"))) invisible(file.remove(file.path(DER, "zcta_la_2010.gpkg")))
st_write(zc, file.path(DER, "zcta_la_2010.gpkg"), quiet = TRUE)
cat("  LA ZCTA polygons written:", nrow(zc), "\n")

hr("9. Percentage uninsured (ACS 2015-2019 5-year, table B27001)")
ad <- file.path(tempdir(), "acs"); dir.create(ad, showWarnings = FALSE)
unzip(file.path(SRC, "20195ca0131000.zip"), exdir = ad)
e <- read.csv(file.path(ad, "e20195ca0131000.txt"), header = FALSE, colClasses = "character")
g <- read.csv(file.path(SRC, "g20195ca.csv"), header = FALSE, colClasses = "character", fileEncoding = "latin1")
# geography file: column 3 = summary level (140 = tract), 5 = logical record number, 49 = GEOID
g <- g[g$V3 == "140" & substr(g$V49, 8, 12) == "06037", ]
a <- merge(data.frame(LOGRECNO = g$V5, Tract = substr(g$V49, 8, 18)), e, by.x = "LOGRECNO", by.y = "V6")
# Table B27001 is the first table in sequence 131 (57 cells from column 7). Within each of the
# nine male and nine female age groups, the third cell is "No health insurance coverage".
cell <- function(k) num(a[[paste0("V", 6 + k)]])
no_coverage <- c(5, 8, 11, 14, 17, 20, 23, 26, 29, 33, 36, 39, 42, 45, 48, 51, 54, 57)
total <- cell(1); uninsured <- Reduce(`+`, lapply(no_coverage, cell))
keep <- which(!is.na(total) & total > 0)
un <- data.frame(Tract = a$Tract[keep], pct_uninsured = round(uninsured[keep] / total[keep] * 100, 4))
check(nrow(un) == 2323 && sum(tr$GEOID %in% un$Tract) == 2322, "2,323 LA tracts with a nonzero denominator; 2,322 of the 2,343 study tracts have a value")
write.csv(un, file.path(DER, "acs_uninsured_2015_2019.csv"), row.names = FALSE)

hr("Done")
print(sessionInfo())
