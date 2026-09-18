# analysis.R: statistical analyses, tables and appendix figures for the article
# "Environmental and Socioeconomic Correlates of Emergency Department Visit Rates for Asthma
#  and Acute Myocardial Infarction in Los Angeles County: A Cross-Sectional Ecological Study"
#
# Inputs : ../data/oct17csv.csv              (the author's merged analysis file; holds the original straight-line exposure)
#          ../data/revisedLAcensustracts.zip (tract geometry, 2010 tracts)
#          ../public-sources/derived/*       (built by 01_prepare_inputs.R and 02_network_distance.R)
# Outputs: output/tables/*.csv, output/spec_log.csv, output/figures/*.png, output/tract_analysis_dataset.csv
# Run from the repository root:
#   Rscript analysis/analysis.R > analysis/output/analysis-log.txt 2>&1
#
# The script stops if any check fails. Every analysis model it fits is printed in the log and
# listed in output/spec_log.csv; auxiliary fits used only for partial R2, variance inflation
# factors and the density-interaction Wald test are not listed separately.
suppressPackageStartupMessages({
  library(sf); library(spdep); library(spatialreg); library(sandwich); library(lmtest)
  library(MASS); library(mgcv); library(ggplot2); library(patchwork); library(ggcorrplot)
})
options(width = 140, stringsAsFactors = FALSE, scipen = 4)
set.seed(20260911)

here <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) dirname(normalizePath(f)) else getwd()
})
PKG <- normalizePath(file.path(here, ".."))
DATA <- file.path(PKG, "data"); DER <- file.path(PKG, "public-sources", "derived")
OUT <- file.path(here, "output"); TAB <- file.path(OUT, "tables"); FIG <- file.path(OUT, "figures")
for (p in c(OUT, TAB, FIG)) dir.create(p, showWarnings = FALSE, recursive = TRUE)

check <- function(ok, msg) { if (!isTRUE(ok)) stop("Check failed: ", msg, call. = FALSE); cat("  ok:", msg, "\n") }
hr   <- function(t) cat("\n", strrep("=", 110), "\n## ", t, "\n", strrep("=", 110), "\n", sep = "")
geoid <- function(x) sprintf("%011.0f", round(as.numeric(x)))
wmean <- function(x, w) sum(w * x) / sum(w)
wsd   <- function(x, w) sqrt(sum(w * (x - wmean(x, w))^2) / sum(w))
SPEC  <- list()                                 # models fitted, written to output/spec_log.csv
log_spec <- function(id, outcome, sample, model, n, term, est, se, lo, hi, p, note = "") {
  num <- function(x) unname(as.numeric(x))
  SPEC[[length(SPEC) + 1]] <<- data.frame(id = id, outcome = outcome, sample = sample, model = model, n = num(n),
                                          term = term, est = num(est), se = num(se), lo = num(lo), hi = num(hi),
                                          p = num(p), note = note, row.names = NULL)
}

hr("0. Data assembly and missing-value correction")
rec <- read.csv(file.path(DATA, "oct17csv.csv"), fileEncoding = "UTF-8-BOM")
rec$GEOID <- geoid(rec$Tract)
c4 <- read.csv(file.path(DER, "ces4_official_la.csv"), colClasses = c(GEOID = "character"))
ex <- read.csv(file.path(DER, "tract_exposures.csv"), colClasses = c(GEOID = "character"))
c3 <- read.csv(file.path(DER, "ces3_ozone_la.csv"), colClasses = c(GEOID = "character"))
un <- read.csv(file.path(DER, "acs_uninsured_2015_2019.csv")); un$GEOID <- geoid(un$Tract)
check(nrow(rec) == 2343 && !anyDuplicated(rec$GEOID), "2,343 tracts in oct17csv.csv")

o <- c4[match(rec$GEOID, c4$GEOID), ]
check(!anyNA(o$GEOID), "every study tract found in the official CalEnviroScreen 4.0 file")
d <- data.frame(GEOID = rec$GEOID, ZIP = as.character(rec$ZIP), pop = rec$Population,
                dist_sl = rec$MEAN_Total_Miles, pm25 = rec$PM2_5, ozone = rec$Ozone,
                asthma_pub = rec$Asthma, ami_pub = rec$Cardiovasc, poverty_pub = rec$Poverty,
                asthma = o$asthma, ami = o$ami, poverty = o$poverty,
                area_sqmi = rec$Shape_Area / 2589988.11,
                child = rec$Child_10, elderly = rec$Elderly_65)
d$density <- ifelse(d$area_sqmi > 0, d$pop / d$area_sqmi, NA)

same <- function(a, b) abs(a - b) < 1e-6
check(all(d$pop == o$population), "population identical to official file")
check(all(same(d$pm25, o$pm25)) && all(same(d$ozone, o$ozone)), "PM2.5 and ozone identical to official file")
for (v in c("asthma", "ami", "poverty")) {
  pub <- d[[paste0(v, "_pub")]]; off <- d[[v]]
  check(all(same(pub[!is.na(off)], off[!is.na(off)])), paste(v, ": every non-missing official value identical"))
  check(all(pub[is.na(off)] == 0), paste(v, ": all", sum(is.na(off)), "officially missing values were carried as 0"))
}
# The two age percentages in the export carry the OEHHA shapefile's four-decimal values; the official
# spreadsheet rounds them to two decimals. They agree except where the export carries 0: the
# 16 tracts without residents (missing in both) and tract 06037930401, which OEHHA's shapefile download
# codes as missing (-999; the packaged tract layer and the export carry 0) although the spreadsheet
# reports its values. The spreadsheet values are used there.
demo <- c(child = "Child_10", elderly = "Elderly_65")
official <- c(child = "child_10", elderly = "elderly_65")
for (v in names(demo)) {
  pub <- rec[[demo[[v]]]]; off <- o[[official[[v]]]]
  differ <- which(!is.na(off) & abs(pub - off) > 0.005 + 1e-9)
  check(all(pub[is.na(off)] == 0) && identical(d$GEOID[differ], "06037930401") && pub[differ] == 0,
        paste0(demo[[v]], ": agrees with the official spreadsheet to two decimals except ", sum(is.na(off)),
               " officially missing values carried as 0 and tract 06037930401 (carried as 0; spreadsheet ", off[differ], ")"))
  d[[v]][is.na(off)] <- NA
  d[[v]][differ] <- off[differ]
}
cat("  zero-population tracts:", sum(d$pop == 0), "| officially missing: asthma", sum(is.na(d$asthma)),
    "| AMI", sum(is.na(d$ami)), "| poverty", sum(is.na(d$poverty)), "\n")
cat("  missing poverty among populated tracts:", sum(d$pop > 0 & is.na(d$poverty)), "\n")
cat("  missing-coded tracts inside the originally submitted <=1.25-mile sample: asthma/AMI",
    sum(d$dist_sl <= 1.25 & is.na(d$asthma)), "| poverty", sum(d$dist_sl <= 1.25 & is.na(d$poverty)), "\n")
check(sum(d$pop == 0) == 16 && sum(is.na(d$asthma)) == 9 && sum(is.na(d$poverty)) == 38, "missing-data counts as expected (16 / 9 / 38)")

m_ex <- match(d$GEOID, ex$GEOID); d <- cbind(d, ex[m_ex, setdiff(names(ex), "GEOID")])
nt <- read.csv(file.path(DER, "tract_network_distance.csv"), colClasses = c(GEOID = "character"))
m_nt <- match(d$GEOID, nt$GEOID)
check(!anyNA(m_nt) && !anyNA(nt$net3_popcentroid[m_nt]), "street-network distances found for every tract")
d$net3_popcentroid <- nt$net3_popcentroid[m_nt]; d$net1_popcentroid <- nt$net1_popcentroid[m_nt]
d$net3_open <- nt$net3_popcentroid_open[m_nt]; d$net3_merged <- nt$net3_popcentroid_merged[m_nt]
d$net3_pna2016 <- nt$net3_popcentroid_pna2016[m_nt]
# Primary exposure: mean street-network distance from the Census 2010 population-weighted tract center to the
# three nearest parks (02_network_distance.R). dist_sl is the original straight-line measure from the tract
# centroid (ArcGIS Find Closest, data/oct17csv.csv), used for the originally submitted models, the 1.25-mile
# comparisons, and the sensitivity rows that name it.
d$dist <- d$net3_popcentroid
d$ozone_ces3 <- c3$ozone_ces3[match(d$GEOID, c3$GEOID)]
d$uninsured  <- un$pct_uninsured[match(d$GEOID, un$GEOID)]
check(all(same(d$recorded_dist3_miles, d$dist_sl)), "exposure table aligned with the original study distances")

d$primary <- d$pop > 0 & complete.cases(d[, c("asthma", "ami", "poverty", "pm25", "ozone", "dist", "dist_sl")])
U <- d[d$primary, ]; rownames(U) <- NULL
cat("  primary sample (residents and complete data):", nrow(U), "tracts\n")
check(nrow(U) == 2305, "primary sample n = 2,305")
U$w <- U$pop
IQR_D <- IQR(U$dist); IQR_SL <- IQR(U$dist_sl)
cat(sprintf("  distance in primary sample: median %.3f, IQR %.3f (Q1 %.3f, Q3 %.3f), SD %.3f, range %.3f-%.3f\n",
            median(U$dist), IQR_D, quantile(U$dist, .25), quantile(U$dist, .75), sd(U$dist), min(U$dist), max(U$dist)))
cat(sprintf("  (street-network distance from the population-weighted center; walking time at 3 mph: median %.1f min, IQR %.1f min)\n",
            median(U$dist) * 20, IQR_D * 20))
cat(sprintf("  original straight-line distance in primary sample: median %.3f, IQR %.3f (Q1 %.3f, Q3 %.3f), SD %.3f, range %.3f-%.3f\n",
            median(U$dist_sl), IQR_SL, quantile(U$dist_sl, .25), quantile(U$dist_sl, .75), sd(U$dist_sl), min(U$dist_sl), max(U$dist_sl)))
cat(sprintf("  Pearson r between the two measures in the primary sample: %.3f\n", cor(U$dist, U$dist_sl)))

hr("1. Check: the originally submitted models are reproduced exactly from the original coding")
pubA <- d[d$dist_sl >= 0 & d$dist_sl <= 1.25, ]
check(nrow(pubA) == 2291, "originally submitted analytic sample n = 2,291")
mp <- lm(asthma_pub ~ dist_sl + pm25 + ozone + poverty_pub, data = pubA)
want <- c(103.349, 13.848, -5.380, -413.834, 0.807)
check(all(abs(round(coef(mp), 3) - want) < 1e-9) && round(summary(mp)$r.squared, 3) == 0.316, "originally submitted Table 4 asthma model (13.848; R2 0.316)")
mpa <- lm(ami_pub ~ dist_sl + pm25 + ozone + poverty_pub, data = pubA)
check(round(coef(mpa)["dist_sl"], 2) == 2.14, "originally submitted AMI model (2.14)")
mpf <- lm(asthma_pub ~ dist_sl + pm25 + ozone + poverty_pub, data = d)
check(round(coef(mpf)["dist_sl"], 2) == 4.93, "originally submitted full-sample sensitivity model (4.93)")
lm_row <- function(m) { cc <- coef(summary(m))["dist_sl", ]; ci <- confint(m)["dist_sl", ]; unname(c(cc[1], cc[2], ci[1], ci[2], cc[4])) }
for (g in list(list("asthma", "originally submitted: <=1.25 mi, missing values coded 0", mp),
               list("ami", "originally submitted: <=1.25 mi, missing values coded 0", mpa),
               list("asthma", "originally submitted: all 2,343 tracts, missing values coded 0", mpf))) {
  r <- lm_row(g[[3]])
  log_spec("G1", g[[1]], g[[2]], "OLS unweighted, conventional SE", nobs(g[[3]]), "original straight-line dist per mile", r[1], r[2], r[3], r[4], r[5])
}
mp10 <- lm(asthma_pub ~ dist_sl + pm25 + ozone + poverty_pub, data = d[d$dist_sl <= 1.0, ])
cat(sprintf("  as originally coded, tracts within 1.0 mile (n = %d): asthma distance coefficient %.3f per mile (1.25-mile threshold: %.3f)\n",
            nobs(mp10), coef(mp10)["dist_sl"], coef(mp)["dist_sl"]))
mpa10 <- lm(ami_pub ~ dist_sl + pm25 + ozone + poverty_pub, data = d[d$dist_sl <= 1.0, ])
cat(sprintf("  as originally coded, tracts within 1.0 mile (n = %d): AMI distance coefficient %.3f per mile (1.25-mile threshold: %.3f)\n",
            nobs(mpa10), coef(mpa10)["dist_sl"], coef(mpa)["dist_sl"]))
log_spec("C9", "asthma", "as originally coded, <=1.0 mi", "OLS unweighted", nobs(mp10), "original straight-line dist per mile", coef(mp10)["dist_sl"], NA, NA, NA, NA)
log_spec("C9", "ami", "as originally coded, <=1.0 mi", "OLS unweighted", nobs(mpa10), "original straight-line dist per mile", coef(mpa10)["dist_sl"], NA, NA, NA, NA)

hr("2. Geometry and spatial weights")
td <- file.path(tempdir(), "tracts"); dir.create(td, showWarnings = FALSE)
unzip(file.path(DATA, "revisedLAcensustracts.zip"), exdir = td)
tr <- st_read(list.files(td, pattern = "[.]shp$", full.names = TRUE)[1], quiet = TRUE)
tr <- st_make_valid(st_transform(tr[, "Tract"], 3310)); tr$GEOID <- geoid(tr$Tract)
geom_of <- function(g) st_geometry(tr)[match(g, tr$GEOID)]
check(all(U$GEOID %in% tr$GEOID), "primary tracts join 1:1 to geometry")

queen_lw <- function(g) {
  nb <- poly2nb(geom_of(g), queen = TRUE)
  list(nb = nb, lw = nb2listw(nb, style = "W", zero.policy = TRUE))
}
knn_lw <- function(dat, k = 6) {
  nb <- knn2nb(knearneigh(cbind(dat$centroid_x, dat$centroid_y), k = k))
  list(nb = nb, lw = nb2listw(nb, style = "W"))
}
QU <- queen_lw(U$GEOID); KU <- knn_lw(U)
cat("  queen weights (primary):", nrow(U), "tracts |", sum(card(QU$nb) == 0), "islands |",
    n.comp.nb(QU$nb)$nc, "sub-graphs | mean links", round(mean(card(QU$nb)), 2), "\n")

f_a <- asthma ~ dist + pm25 + ozone + poverty
f_m <- ami ~ dist + pm25 + ozone + poverty
m0 <- lm(f_a, data = U, weights = w)             # the Table 4 asthma model, logged in section 6
# Moran's I for regression residuals as in spdep::lm.morantest: weighted residuals sqrt(w) * e, not re-centered
e <- weighted.residuals(m0); W <- as(QU$lw, "CsparseMatrix"); keep <- card(QU$nb) > 0
I_hand <- (sum(keep) / sum(W)) * as.numeric(t(e) %*% (W %*% e)) / sum(e^2)
I_pkg <- lm.morantest(m0, QU$lw, zero.policy = TRUE)$estimate[1]
cat(sprintf("  Moran's I of weighted-OLS asthma residuals (the value reported in section 6): hand %.4f | spdep::lm.morantest %.4f\n", I_hand, I_pkg))
check(abs(I_hand - I_pkg) < 0.005, "hand-computed Moran's I equals package value")
I_scr <- moran.test(sample(e), QU$lw, zero.policy = TRUE)$estimate[1]
cat(sprintf("  negative control, scrambled residuals: I = %.4f\n", I_scr))
check(abs(I_scr) < 0.05, "scrambled residuals show no autocorrelation")

# Model helpers
zip_vc <- function(m, dat) vcovCL(m, cluster = dat$ZIP, type = "HC1")
term_row <- function(m, term, V, df, scale = 1) {
  est <- unname(coef(m)[term]); se <- unname(sqrt(V[term, term])); t <- est / se
  p <- if (is.finite(df)) 2 * pt(-abs(t), df) else 2 * pnorm(-abs(t))
  q <- if (is.finite(df)) qt(0.975, df) else qnorm(0.975)
  c(est = scale * est, se = scale * se, lo = scale * (est - q * se), hi = scale * (est + q * se), p = p)
}
partial_r2 <- function(m, dat, term) {
  red <- update(m, as.formula(paste(". ~ . -", term)), data = dat)
  (deviance(red) - deviance(m)) / deviance(red)
}
fit_ols <- function(f, dat, weighted = TRUE, term = "dist", scale = 1) {
  dat$w <- if (weighted) dat$pop else rep(1, nrow(dat))
  m <- lm(f, data = dat, weights = w)
  G <- length(unique(dat$ZIP))
  zc <- term_row(m, term, zip_vc(m, dat), G - 1, scale)
  h1 <- term_row(m, term, vcovHC(m, type = "HC1"), m$df.residual, scale)
  y <- model.response(model.frame(m)); x <- dat[[term]]
  list(m = m, n = nobs(m), zip = zc, hc1 = h1, r2 = summary(m)$r.squared, adj = summary(m)$adj.r.squared,
       pr2 = partial_r2(m, dat, term),
       std = unname(coef(m)[term]) * wsd(x, dat$w) / wsd(y, dat$w))
}
fit_sem <- function(f, dat, lwobj, weighted = TRUE, term = "dist", scale = 1) {
  dat$w <- if (weighted) dat$pop else rep(1, nrow(dat))
  m <- errorsarlm(f, data = dat, listw = lwobj$lw, weights = w, zero.policy = TRUE)
  s <- summary(m)
  est <- unname(coef(m)[term]); se <- unname(s$Coef[term, "Std. Error"])
  lr <- tryCatch(unname(LR1.Sarlm(m)$p.value[1]), error = function(e) NA)
  mi <- moran.test(residuals(m), lwobj$lw, zero.policy = TRUE)
  list(m = m, n = length(residuals(m)),
       row = c(est = scale * est, se = scale * se, lo = scale * (est - 1.96 * se), hi = scale * (est + 1.96 * se),
               p = 2 * pnorm(-abs(est / se))),
       lambda = m$lambda, lr_p = lr, resid_I = unname(mi$estimate[1]), resid_I_p = mi$p.value)
}
fmt <- function(v, k = 3) formatC(v, digits = k, format = "f")
pfmt <- function(p) ifelse(p < 0.001, formatC(p, format = "e", digits = 1), formatC(p, format = "f", digits = 3))

hr("3. Table 1: descriptive statistics, primary sample")
vars <- c(dist = "Mean street-network distance to three nearest parks (miles)",
          walk3 = "Walking time to three nearest parks at 3 mph (minutes)",
          dist_sl = "Mean straight-line distance to three nearest parks, from tract centroid (miles)", asthma = "Asthma ED visit rate (per 10,000)",
          ami = "AMI ED visit rate (per 10,000)", pm25 = "PM2.5 (ug/m3)", ozone = "Ozone (ppm)",
          poverty = "Residents below twice the federal poverty level (%)", pop = "Tract population",
          density = "Population density (per sq mi)")
U$walk3 <- U$dist * 20
t1 <- do.call(rbind, lapply(names(vars), function(v) {
  x <- U[[v]]
  data.frame(variable = vars[[v]], n = length(x), mean = mean(x), sd = sd(x), median = median(x),
             q1 = unname(quantile(x, .25)), q3 = unname(quantile(x, .75)), min = min(x), max = max(x),
             pop_weighted_mean = wmean(x, U$pop))
}))
print(t1, digits = 4); write.csv(t1, file.path(TAB, "table1_descriptives.csv"), row.names = FALSE)

hr("4. Table 2: tracts within vs beyond 1.25 miles by the original straight-line measure (primary sample)")
inn <- U$dist_sl <= 1.25
cat("  within 1.25 mi:", sum(inn), "| beyond:", sum(!inn), "\n")
t2 <- do.call(rbind, lapply(setdiff(names(vars), "walk3"), function(v) {
  a <- U[[v]][inn]; b <- U[[v]][!inn]
  data.frame(variable = vars[[v]], n_within = length(a), mean_within = mean(a), sd_within = sd(a),
             median_within = median(a), q1_within = unname(quantile(a, .25)), q3_within = unname(quantile(a, .75)),
             n_beyond = length(b), mean_beyond = mean(b), sd_beyond = sd(b),
             median_beyond = median(b), q1_beyond = unname(quantile(b, .25)), q3_beyond = unname(quantile(b, .75)),
             smd = (mean(a) - mean(b)) / sqrt((var(a) + var(b)) / 2),
             wilcoxon_p = suppressWarnings(wilcox.test(a, b, exact = FALSE)$p.value))
}))
print(t2, digits = 4); write.csv(t2, file.path(TAB, "table2_within_vs_beyond_1p25mi.csv"), row.names = FALSE)
cat("  original definition (all 2,343; zeros as coded): excluded tracts =", sum(d$dist_sl > 1.25), "\n")

hr("5. Table 3: bivariate associations of park distance (population-weighted, ZIP-clustered)")
t3 <- do.call(rbind, lapply(c("asthma", "ami", "pm25", "ozone", "poverty"), function(y) {
  r <- fit_ols(as.formula(paste(y, "~ dist")), U)
  log_spec("T3", y, "primary", "WLS bivariate", r$n, "dist per mile", r$zip["est"], r$zip["se"], r$zip["lo"], r$zip["hi"], r$zip["p"])
  data.frame(outcome = y, n = r$n, per_mile = r$zip["est"], lo = r$zip["lo"], hi = r$zip["hi"], p_zip = r$zip["p"],
             p_hc1 = r$hc1["p"], per_iqr = r$zip["est"] * IQR_D, r2 = r$r2)
}))
t3u <- do.call(rbind, lapply(c("asthma", "ami", "pm25", "ozone", "poverty"), function(y) {
  r <- fit_ols(as.formula(paste(y, "~ dist")), U, weighted = FALSE)
  log_spec("T3u", y, "primary", "OLS bivariate unweighted (ZIP-clustered)", r$n, "dist per mile", r$zip["est"], r$zip["se"], r$zip["lo"], r$zip["hi"], r$zip["p"])
  data.frame(outcome = y, unweighted_per_mile = r$zip["est"], unweighted_p_zip = r$zip["p"], unweighted_r2 = r$r2)
}))
t3 <- cbind(t3, t3u[, -1]); rownames(t3) <- NULL
print(t3, digits = 4); write.csv(t3, file.path(TAB, "table3_bivariate.csv"), row.names = FALSE)
cat("\n  The same weighted bivariate models with the original straight-line measure (not tabulated):\n")
for (y in c("asthma", "ami", "pm25", "ozone", "poverty")) {
  r <- fit_ols(as.formula(paste(y, "~ dist_sl")), U, term = "dist_sl")
  cat(sprintf("    %-8s %.4f per mile (95%% CI %.4f to %.4f), ZIP-cl p %s, R2 %.4f\n", y, r$zip["est"], r$zip["lo"], r$zip["hi"], pfmt(r$zip["p"]), r$r2))
  log_spec("T3sl", y, "primary", "WLS bivariate, original straight-line measure", r$n, "original straight-line dist per mile", r$zip["est"], r$zip["se"], r$zip["lo"], r$zip["hi"], r$zip["p"])
}
cat("\n  Distance-poverty association by sample definition (three definitions), original straight-line measure, unweighted OLS, conventional SEs:\n")
for (sp in list(list("all 2,343 tracts, as originally coded (missing poverty = 0)", d, "poverty_pub"),
                list("tracts with a reported poverty value", d[!is.na(d$poverty), ], "poverty"),
                list("as originally coded, tracts within 1.25 miles", d[d$dist_sl <= 1.25, ], "poverty_pub"))) {
  mm <- lm(as.formula(paste(sp[[3]], "~ dist_sl")), data = sp[[2]])
  cc <- coef(summary(mm))["dist_sl", ]; ci95 <- confint(mm)["dist_sl", ]
  cat(sprintf("    %-58s n %4d: %.3f per mile (95%% CI %.3f to %.3f), p %s\n", sp[[1]], nobs(mm), cc[1], ci95[1], ci95[2], pfmt(cc[4])))
  log_spec("C8", "poverty", sp[[1]], "OLS bivariate, conventional SE", nobs(mm), "original straight-line dist per mile", cc[1], cc[2], ci95[1], ci95[2], cc[4])
}

hr("6. Table 4: main models (primary sample)")
main <- list(); t4 <- list()
for (oc in c("asthma", "ami")) {
  f <- if (oc == "asthma") f_a else f_m
  A  <- fit_ols(f, U)                      # weighted OLS
  B  <- fit_sem(f, U, QU)                  # weighted SEM, queen
  Au <- fit_ols(f, U, weighted = FALSE)
  Bu <- fit_sem(f, U, QU, weighted = FALSE)
  mi <- lm.morantest(A$m, QU$lw, zero.policy = TRUE)
  ha <- tryCatch(Hausman.test(B$m), error = function(e) NULL)
  cat(sprintf("\n[%s] weighted OLS: dist per mile %.3f (ZIP-cl 95%% CI %.3f to %.3f, p %s; HC1 p %s) | per IQR %.3f | R2 %.3f | partial R2 %.4f | std beta %.3f | resid Moran's I %.3f (p %s)\n",
              oc, A$zip["est"], A$zip["lo"], A$zip["hi"], pfmt(A$zip["p"]), pfmt(A$hc1["p"]), A$zip["est"] * IQR_D,
              A$r2, A$pr2, A$std, mi$estimate[1], pfmt(mi$p.value)))
  cat(sprintf("[%s] weighted SEM: dist per mile %.3f (95%% CI %.3f to %.3f, p %s) | per IQR %.3f | lambda %.3f (LR p %s) | resid Moran's I %.3f (p %s) | Hausman p %s\n",
              oc, B$row["est"], B$row["lo"], B$row["hi"], pfmt(B$row["p"]), B$row["est"] * IQR_D, B$lambda,
              pfmt(B$lr_p), B$resid_I, pfmt(B$resid_I_p), if (is.null(ha)) "NA" else pfmt(ha$p.value)))
  cat(sprintf("[%s] unweighted OLS: %.3f (ZIP-cl p %s) | unweighted SEM: %.3f (p %s, lambda %.3f)\n",
              oc, Au$zip["est"], pfmt(Au$zip["p"]), Bu$row["est"], pfmt(Bu$row["p"]), Bu$lambda))
  # full coefficient table, both primary models
  G <- length(unique(U$ZIP)); VA <- zip_vc(A$m, U); sB <- summary(B$m)
  for (tm in c("(Intercept)", "dist", "pm25", "ozone", "poverty")) {
    sc <- if (tm == "ozone") 0.01 else 1
    ra <- term_row(A$m, tm, VA, G - 1, sc)
    eb <- unname(coef(B$m)[tm]); sb <- unname(sB$Coef[tm, "Std. Error"])
    rb <- c(est = sc * eb, se = sc * sb, lo = sc * (eb - 1.96 * sb), hi = sc * (eb + 1.96 * sb), p = 2 * pnorm(-abs(eb / sb)))
    t4[[length(t4) + 1]] <- data.frame(outcome = oc, term = tm, unit = ifelse(tm == "ozone", "per 0.01 ppm", ifelse(tm == "dist", "per mile", "per unit")),
                                       ols_est = ra["est"], ols_lo = ra["lo"], ols_hi = ra["hi"], ols_p = ra["p"],
                                       sem_est = rb["est"], sem_lo = rb["lo"], sem_hi = rb["hi"], sem_p = rb["p"])
    log_spec("T4", oc, "primary", "WLS (ZIP-clustered)", A$n, tm, ra["est"], ra["se"], ra["lo"], ra["hi"], ra["p"], ifelse(tm == "ozone", "per 0.01 ppm", ""))
    log_spec("T4", oc, "primary", "weighted SEM queen", B$n, tm, rb["est"], rb["se"], rb["lo"], rb["hi"], rb["p"], ifelse(tm == "ozone", "per 0.01 ppm", ""))
  }
  log_spec("T4", oc, "primary", "OLS unweighted (ZIP-clustered)", Au$n, "dist", Au$zip["est"], Au$zip["se"], Au$zip["lo"], Au$zip["hi"], Au$zip["p"])
  log_spec("T4", oc, "primary", "SEM unweighted queen", Bu$n, "dist", Bu$row["est"], Bu$row["se"], Bu$row["lo"], Bu$row["hi"], Bu$row["p"])
  main[[oc]] <- list(A = A, B = B, Au = Au, Bu = Bu, moran = mi, hausman = ha)
}
cat("\n  Standardized coefficients (weighted SDs) and partial R2 of every predictor, weighted OLS:\n")
for (oc in c("asthma", "ami")) {
  A <- main[[oc]]$A$m; y <- U[[oc]]
  sb <- sapply(c("dist", "pm25", "ozone", "poverty"), function(v) unname(coef(A)[v]) * wsd(U[[v]], U$pop) / wsd(y, U$pop))
  pr <- sapply(c("dist", "pm25", "ozone", "poverty"), function(v) partial_r2(A, U, v))
  cat(sprintf("  [%s] std beta: %s | partial R2: %s\n", oc, paste(names(sb), round(sb, 3), collapse = " "),
              paste(names(pr), round(pr, 4), collapse = " ")))
  for (v in names(sb)) log_spec("T4std", oc, "primary", "WLS standardized coefficient", nobs(A), v, sb[v], NA, NA, NA, NA, paste("partial R2", round(pr[v], 4)))
}
t4 <- do.call(rbind, t4); rownames(t4) <- NULL
print(t4, digits = 4); write.csv(t4, file.path(TAB, "table4a_main_models.csv"), row.names = FALSE)
t4s <- do.call(rbind, lapply(names(main), function(oc) with(main[[oc]], data.frame(
  outcome = oc, n = A$n, iqr_miles = IQR_D,
  ols_per_iqr = A$zip["est"] * IQR_D, ols_lo = A$zip["lo"] * IQR_D, ols_hi = A$zip["hi"] * IQR_D, ols_p_zip = A$zip["p"], ols_p_hc1 = A$hc1["p"],
  ols_r2 = A$r2, ols_partial_r2 = A$pr2, ols_std_beta = A$std, ols_resid_moran = moran$estimate[1], ols_resid_moran_p = moran$p.value,
  sem_per_iqr = B$row["est"] * IQR_D, sem_lo = B$row["lo"] * IQR_D, sem_hi = B$row["hi"] * IQR_D, sem_p = B$row["p"],
  sem_lambda = B$lambda, sem_lr_p = B$lr_p, sem_resid_moran = B$resid_I,
  hausman_p = if (is.null(hausman)) NA else hausman$p.value,
  unw_ols_per_iqr = Au$zip["est"] * IQR_D, unw_ols_p = Au$zip["p"], unw_sem_per_iqr = Bu$row["est"] * IQR_D, unw_sem_p = Bu$row["p"]))))
rownames(t4s) <- NULL; print(t4s, digits = 4); write.csv(t4s, file.path(TAB, "table4a_summary.csv"), row.names = FALSE)

hr("7. Table 5: within vs between ZIP codes")
vshare <- function(x, g, w = NULL) {
  if (is.null(w)) w <- rep(1, length(x))
  gm <- ave(x * w, g, FUN = sum) / ave(w, g, FUN = sum)
  sum(w * (x - gm)^2) / sum(w * (x - wmean(x, w))^2)
}
t4b1 <- data.frame(variable = c("asthma", "ami", "dist", "pm25", "ozone", "poverty"))
t4b1$within_zip_share_unweighted <- sapply(t4b1$variable, function(v) vshare(U[[v]], U$ZIP))
t4b1$within_zip_share_weighted   <- sapply(t4b1$variable, function(v) vshare(U[[v]], U$ZIP, U$pop))
print(t4b1, digits = 3)
cat("  ZIP codes:", length(unique(U$ZIP)), "| tracts per ZIP: median", median(table(U$ZIP)), "\n")
for (v in c("dist", "pm25", "ozone", "poverty")) {
  U[[paste0(v, "_b")]] <- ave(U[[v]] * U$pop, U$ZIP, FUN = sum) / ave(U$pop, U$ZIP, FUN = sum)
  U[[paste0(v, "_w")]] <- U[[v]] - U[[paste0(v, "_b")]]
}
t4b2 <- list()
for (oc in c("asthma", "ami")) {
  f <- as.formula(paste(oc, "~ dist_w + dist_b + pm25_w + pm25_b + ozone_w + ozone_b + poverty_w + poverty_b"))
  m <- lm(f, data = U, weights = pop); V <- zip_vc(m, U); G <- length(unique(U$ZIP))
  for (tm in c("dist_w", "dist_b", "pm25_w", "pm25_b")) {
    r <- term_row(m, tm, V, G - 1)
    t4b2[[length(t4b2) + 1]] <- data.frame(outcome = oc, term = tm, est_per_unit = r["est"], lo = r["lo"], hi = r["hi"], p = r["p"],
                                           est_per_iqr = if (grepl("^dist", tm)) r["est"] * IQR_D else NA)
    log_spec("T4B", oc, "primary", "WLS within/between ZIP", nobs(m), tm, r["est"], r["se"], r["lo"], r["hi"], r["p"])
  }
  mf <- lm(update(f_a, paste(oc, "~ . + factor(ZIP)")), data = U, weights = pop)
  r <- term_row(mf, "dist", zip_vc(mf, U), G - 1)
  cat(sprintf("  [%s] ZIP fixed effects: dist per mile %.3f (95%% CI %.3f to %.3f, p %s)\n", oc, r["est"], r["lo"], r["hi"], pfmt(r["p"])))
  log_spec("T4B", oc, "primary", "WLS ZIP fixed effects", nobs(mf), "dist", r["est"], r["se"], r["lo"], r["hi"], r["p"])
}
t4b2 <- do.call(rbind, t4b2); rownames(t4b2) <- NULL; print(t4b2, digits = 4)
write.csv(t4b1, file.path(TAB, "table4b_variance_shares.csv"), row.names = FALSE)
write.csv(t4b2, file.path(TAB, "table4b_within_between.csv"), row.names = FALSE)

hr("8. Tables 6 and 7: sensitivity analyses")
S <- list()
add_sens <- function(id, label, dat, xvar = "dist", xscale = NULL, weighted = TRUE, sem = TRUE, lwobj = NULL,
                     extra = "", ozvar = "ozone", note = "") {
  dat$x <- dat[[xvar]]; dat$oz <- dat[[ozvar]]
  if (is.null(xscale)) xscale <- IQR(dat$x)
  for (oc in c("asthma", "ami")) {
    f <- as.formula(paste(oc, "~ x + pm25 + oz + poverty", extra))
    r <- fit_ols(f, dat, weighted = weighted, term = "x", scale = xscale)
    pv <- term_row(r$m, "poverty", zip_vc(r$m, dat), length(unique(dat$ZIP)) - 1)
    row <- data.frame(id, label, outcome = oc, n = r$n, scale = xscale,
                      ols_est = r$zip["est"], ols_lo = r$zip["lo"], ols_hi = r$zip["hi"], ols_p = r$zip["p"],
                      r2 = r$r2, partial_r2 = r$pr2, std_beta = r$std,
                      sem_est = NA, sem_lo = NA, sem_hi = NA, sem_p = NA, lambda = NA,
                      pov_ols = pv["est"], pov_ols_p = pv["p"], pov_sem = NA, pov_sem_p = NA, note, row.names = NULL)
    log_spec(id, oc, label, ifelse(weighted, "WLS (ZIP-clustered)", "OLS (ZIP-clustered)"), r$n, paste("exposure per", signif(xscale, 3)),
             r$zip["est"], r$zip["se"], r$zip["lo"], r$zip["hi"], r$zip["p"], note)
    if (sem) {
      lwo <- if (is.null(lwobj)) queen_lw(dat$GEOID) else lwobj
      s <- tryCatch(fit_sem(f, dat, lwo, weighted = weighted, term = "x", scale = xscale), error = function(e) NULL)
      if (!is.null(s)) {
        row[, c("sem_est", "sem_lo", "sem_hi", "sem_p", "lambda")] <- c(s$row[c("est", "lo", "hi", "p")], s$lambda)
        row[, c("pov_sem", "pov_sem_p")] <- summary(s$m)$Coef["poverty", c(1, 4)]
        log_spec(id, oc, label, ifelse(weighted, "weighted SEM", "SEM unweighted"), s$n, paste("exposure per", signif(xscale, 3)),
                 s$row["est"], s$row["se"], s$row["lo"], s$row["hi"], s$row["p"], note)
      }
    }
    S[[length(S) + 1]] <<- row
  }
}
# S0 is the primary analysis (street-network distance). S24 and the originally submitted rows (S1, S2, S23) use the original
# straight-line measure (dist_sl) with its own IQR; the threshold rows (S4, S5) restrict the sample by the original measure.
add_sens("S0", "Primary analysis (population-weighted)", U, xscale = IQR_D, lwobj = QU)
add_sens("S24", "Straight-line measure: distance from geographic centroid to three nearest parks", U, xvar = "dist_sl", xscale = IQR_SL, lwobj = QU)
pubF <- d; pubF$asthma <- d$asthma_pub; pubF$ami <- d$ami_pub; pubF$poverty <- d$poverty_pub
add_sens("S1", "As originally submitted: straight-line measure, all 2,343 tracts, missing values coded 0, unweighted", pubF, xvar = "dist_sl", xscale = IQR_SL, weighted = FALSE)
add_sens("S2", "As originally submitted: straight-line measure, <=1.25 mi (n = 2,291), missing values coded 0, unweighted", pubF[pubF$dist_sl <= 1.25, ], xvar = "dist_sl", xscale = IQR_SL, weighted = FALSE)
add_sens("S23", "As originally submitted: straight-line measure, <=1.0 mi (n = 2,237), missing values coded 0, unweighted", pubF[pubF$dist_sl <= 1.0, ], xvar = "dist_sl", xscale = IQR_SL, weighted = FALSE)
add_sens("S3", "Primary sample, unweighted", U, xscale = IQR_D, weighted = FALSE, lwobj = QU)
add_sens("S4", "Primary sample restricted to tracts within 1.25 mi by the straight-line measure", U[U$dist_sl <= 1.25, ], xscale = IQR_D)
add_sens("S5", "Primary sample restricted to tracts within 1.0 mi by the straight-line measure", U[U$dist_sl <= 1.0, ], xscale = IQR_D)
add_sens("S6", "Density >= 1,000 persons per sq mi", U[which(U$density >= 1000), ], xscale = IQR_D)
tert <- cut(U$density, quantile(U$density, c(0, 1/3, 2/3, 1)), include.lowest = TRUE, labels = c("T1", "T2", "T3"))
for (tt in levels(tert)) add_sens(paste0("S7", tt), paste("Density tertile", tt, "(post hoc)"), U[which(tert == tt), ], xscale = IQR_D, sem = FALSE)
for (oc in c("asthma", "ami")) {
  U$tert <- tert
  mi1 <- lm(as.formula(paste(oc, "~ dist * tert + pm25 + ozone + poverty")), data = U, weights = pop)
  mi0 <- lm(as.formula(paste(oc, "~ dist + tert + pm25 + ozone + poverty")), data = U, weights = pop)
  wt <- waldtest(mi0, mi1, vcov = zip_vc(mi1, U))
  cat(sprintf("  [%s] distance x density-tertile interaction: joint Wald p (ZIP-clustered) = %s\n", oc, pfmt(wt[2, "Pr(>F)"])))
  log_spec("S7int", oc, "primary", "WLS dist x density tertile", nobs(mi1), "joint interaction", NA, NA, NA, NA, wt[2, "Pr(>F)"])
}
add_sens("S8", "Six-nearest-neighbor spatial weights", U, xscale = IQR_D, lwobj = KU)
av <- substr(U$ZIP, 1, 3) == "935"
cat("  Antelope Valley (ZIP 935xx) tracts in primary sample:", sum(av), "\n")
add_sens("S9", "Excluding Antelope Valley ZIP codes 935xx (post hoc diagnostic)", U[!av, ], xscale = IQR_D)
U$log2dist <- log2(U$dist)
add_sens("S10", "log2(distance), per doubling", U, xvar = "log2dist", xscale = 1, lwobj = QU)
best <- if (cor(U$dist3_geodesic, U$dist_sl) >= cor(U$dist3_planar, U$dist_sl)) "dist3_geodesic" else "dist3_planar"
cat("  recomputed straight-line distance used for comparison with the straight-line measure:", best, "\n")
add_sens("S11", "Straight-line distance from geographic centroid recomputed from geometry (same park layer)", U, xvar = best, lwobj = QU)
add_sens("S12", "Straight-line distance from population-weighted centroid (Census 2010)", U, xvar = "dist3_popcentroid", lwobj = QU)
add_sens("S22", "Street-network distance from population-weighted centroid to the single nearest park", U, xvar = "net1_popcentroid", lwobj = QU)
add_sens("S25", "Street-network distance to three nearest open-access parks (Restricted and Unknown access excluded)", U, xvar = "net3_open", lwobj = QU)
add_sens("S26", "Street-network distance to three nearest parks with overlapping polygons merged", U, xvar = "net3_merged", lwobj = QU)
add_sens("S13", "Straight-line distance from geographic centroid to the single nearest park", U, xvar = "dist1_planar", lwobj = QU)
add_sens("S14", "Park acres within 0.5 mile (more acres = better access)", U, xvar = "acres_0p5mi", lwobj = QU)
add_sens("S15", "Park acres within 1 mile (more acres = better access)", U, xvar = "acres_1mi", lwobj = QU)
add_sens("S16", "Straight-line distance from geographic centroid to the 2016 county park inventory (contemporaneous with outcomes)", U, xvar = "dist3_pna2016", lwobj = QU)
add_sens("S27", "Street-network distance from population-weighted centroid to three nearest parks in the 2016 county park inventory (parks, roads and outcomes in one period)", U, xvar = "net3_pna2016", lwobj = QU)
add_sens("S17", "Share of tract within a 2016 half-mile walkshed (pedestrian network; more = better access)", U, xvar = "walkshed2016_share", lwobj = QU)
add_sens("S18", "Ozone from CalEnviroScreen 3.0 (2012-2014)", U, xscale = IQR_D, ozvar = "ozone_ces3", lwobj = QU)
U$ozone_mean34 <- (U$ozone + U$ozone_ces3) / 2
add_sens("S19", "Ozone = mean of CalEnviroScreen 3.0 and 4.0 (brackets 2015-2017)", U, xscale = IQR_D, ozvar = "ozone_mean34", lwobj = QU)
T5 <- do.call(rbind, S); rownames(T5) <- NULL
print(T5[, c("id", "outcome", "n", "scale", "ols_est", "ols_lo", "ols_hi", "ols_p", "r2", "std_beta", "sem_est", "sem_lo", "sem_hi", "sem_p", "lambda")], digits = 3)
write.csv(T5, file.path(TAB, "table5_sensitivity.csv"), row.names = FALSE)
cat(sprintf("\n  Poverty across all sensitivity models: OLS coefficient min %.3f (largest p %s); SEM coefficient min %.3f (largest p %s)\n",
            min(T5$pov_ols), pfmt(max(T5$pov_ols_p)), min(T5$pov_sem, na.rm = TRUE), pfmt(max(T5$pov_sem_p, na.rm = TRUE))))

cat("\n  Exposure-metric agreement with the primary (street-network) metric (Pearson r, primary sample):\n")
print(round(cor(U[, c("dist", "dist_sl", best, "dist3_popcentroid", "net1_popcentroid", "net3_open", "net3_merged", "net3_pna2016", "dist1_planar", "acres_0p5mi", "acres_1mi", "dist3_pna2016", "walkshed2016_share")])[1, ], 3))
cat("  Agreement with the original straight-line metric (Pearson r, primary sample):\n")
print(round(cor(U[, c("dist_sl", "dist", best, "dist3_popcentroid", "dist1_planar", "acres_0p5mi", "acres_1mi", "dist3_pna2016", "walkshed2016_share")])[1, ], 3))
cat(sprintf("  street-network distance (three nearest parks) vs straight-line distance from the same population-weighted centers: r %.3f; median ratio %.2f; walking time at 3 mph: median %.1f min, IQR %.1f min\n",
            cor(U$net3_popcentroid, U$dist3_popcentroid), median(U$net3_popcentroid / U$dist3_popcentroid),
            median(U$net3_popcentroid) * 20, IQR(U$net3_popcentroid) * 20))
cat(sprintf("  CES 3.0 vs CES 4.0 ozone: Pearson r %.3f, Spearman r %.3f, mean difference (4.0 - 3.0) %.4f ppm\n",
            cor(U$ozone, U$ozone_ces3), cor(U$ozone, U$ozone_ces3, method = "spearman"), mean(U$ozone - U$ozone_ces3)))

hr("9. Table 8: ZIP-level asthma model (unsmoothed counts at the level where visits are recorded)")
rel <- read.csv(file.path(DER, "zcta_tract_rel_la.csv"), colClasses = c(ZCTA5 = "character", GEOID = "character"))
zip <- read.csv(file.path(DER, "cdph_asthma_zip_la.csv"), colClasses = c(zip = "character"))
rr <- merge(rel, U[, c("GEOID", "dist", "pm25", "ozone", "poverty", "child", "elderly")], by = "GEOID")
rr <- rr[rr$POPPT > 0, ]
agg <- do.call(rbind, lapply(split(rr, rr$ZCTA5), function(z) {
  w <- z$POPPT
  data.frame(ZCTA5 = z$ZCTA5[1], ZPOP = z$ZPOP[1], covered = sum(w) / z$ZPOP[1],
             dist = wmean(z$dist, w), pm25 = wmean(z$pm25, w), ozone = wmean(z$ozone, w),
             poverty = wmean(z$poverty, w), child = wmean(z$child, w), elderly = wmean(z$elderly, w))
}))
Z <- merge(agg, zip[zip$years == 3, ], by.x = "ZCTA5", by.y = "zip")
Z <- Z[Z$covered >= 0.9, ]
cat("  ZIP/ZCTAs analyzed (3 years of counts, >=90% of population in primary tracts):", nrow(Z),
    "| total visits 2015-2017:", sum(Z$visits_2015_2017), "\n")
IQR_Z <- IQR(Z$dist)
nbm <- glm.nb(visits_2015_2017 ~ dist + pm25 + ozone + poverty + child + elderly + offset(log(3 * ZPOP)), data = Z)
b <- coef(nbm)["dist"]; se <- sqrt(vcov(nbm)["dist", "dist"])
cat(sprintf("  negative binomial: rate ratio per ZCTA IQR of distance (%.3f mi) = %.3f (95%% CI %.3f to %.3f), p %s | theta %.2f\n",
            IQR_Z, exp(b * IQR_Z), exp((b - 1.96 * se) * IQR_Z), exp((b + 1.96 * se) * IQR_Z), pfmt(2 * pnorm(-abs(b / se))), nbm$theta))
log_spec("Z1", "asthma visits", "ZCTA", "negative binomial + offset", nrow(Z), "log RR per mile", b, se, b - 1.96 * se, b + 1.96 * se, 2 * pnorm(-abs(b / se)))
zp <- coef(summary(nbm))["poverty", ]
cat(sprintf("  negative binomial, poverty: log rate ratio %.4f per percentage point, p %s\n", zp[1], pfmt(zp[4])))
zc <- st_read(file.path(DER, "zcta_la_2010.gpkg"), quiet = TRUE)
zc <- zc[match(Z$ZCTA5, zc$ZCTA5), ]
check(!anyNA(zc$ZCTA5) && identical(zc$ZCTA5, Z$ZCTA5), "every analyzed ZCTA has a polygon, in the same order")
znb <- poly2nb(zc, queen = TRUE)
if (any(card(znb) == 0)) {                     # link any island to its nearest ZCTA
  cc <- st_coordinates(st_centroid(st_geometry(zc)))
  nn <- knn2nb(knearneigh(cc, k = 1))
  for (i in which(card(znb) == 0)) { j <- nn[[i]]; znb[[i]] <- as.integer(j); znb[[j]] <- sort(unique(c(znb[[j]][znb[[j]] > 0], i))) }
}
zlw <- nb2listw(znb, style = "W")
zmi <- moran.test(residuals(nbm, type = "pearson"), zlw)
cat(sprintf("  Moran's I of Pearson residuals (ZCTA contiguity): %.3f, p %s\n", zmi$estimate[1], pfmt(zmi$p.value)))
nbl <- lapply(seq_along(znb), function(i) Z$ZCTA5[znb[[i]]]); names(nbl) <- Z$ZCTA5
Z$fz <- factor(Z$ZCTA5, levels = names(nbl))
gm <- gam(visits_2015_2017 ~ dist + pm25 + ozone + poverty + child + elderly + s(fz, bs = "mrf", xt = list(nb = nbl)) +
            offset(log(3 * ZPOP)), family = nb(), data = Z, method = "REML")
gb <- coef(gm)["dist"]; gse <- sqrt(vcov(gm)["dist", "dist"])
cat(sprintf("  negative binomial + ZCTA Markov random field: rate ratio per IQR = %.3f (95%% CI %.3f to %.3f), p %s\n",
            exp(gb * IQR_Z), exp((gb - 1.96 * gse) * IQR_Z), exp((gb + 1.96 * gse) * IQR_Z), pfmt(2 * pnorm(-abs(gb / gse)))))
log_spec("Z2", "asthma visits", "ZCTA", "negative binomial + MRF spatial effect", nrow(Z), "log RR per mile", gb, gse, gb - 1.96 * gse, gb + 1.96 * gse, 2 * pnorm(-abs(gb / gse)))
gp <- summary(gm)$p.table["poverty", ]
cat(sprintf("  negative binomial + MRF, poverty: log rate ratio %.4f per percentage point, p %s\n", gp[1], pfmt(gp[4])))
zl <- lm(mean_age_adj_rate ~ dist + pm25 + ozone + poverty, data = Z, weights = ZPOP)
zr <- coeftest(zl, vcov = vcovHC(zl, type = "HC1"))["dist", ]
qz <- qt(0.975, zl$df.residual)               # t-based CI, matching the t-based HC1 p-value (as for the tract-level HC1 rows)
cat(sprintf("  weighted linear model of the unsmoothed age-adjusted ZIP rate: %.2f per mile (%.2f per IQR), HC1 p %s\n",
            zr[1], zr[1] * IQR_Z, pfmt(zr[4])))
log_spec("Z3", "asthma age-adjusted rate", "ZCTA", "WLS (HC1)", nobs(zl), "dist per mile", zr[1], zr[2], zr[1] - qz * zr[2], zr[1] + qz * zr[2], zr[4])
tz <- data.frame(model = c("NB + offset", "NB + offset + MRF", "WLS on unsmoothed age-adjusted rate"),
                 n_zcta = nrow(Z), iqr_miles = IQR_Z,
                 estimate = c(exp(b * IQR_Z), exp(gb * IQR_Z), zr[1] * IQR_Z),
                 lo = c(exp((b - 1.96 * se) * IQR_Z), exp((gb - 1.96 * gse) * IQR_Z), (zr[1] - qz * zr[2]) * IQR_Z),
                 hi = c(exp((b + 1.96 * se) * IQR_Z), exp((gb + 1.96 * gse) * IQR_Z), (zr[1] + qz * zr[2]) * IQR_Z),
                 p = c(2 * pnorm(-abs(b / se)), 2 * pnorm(-abs(gb / gse)), zr[4]),
                 scale = c("rate ratio per IQR", "rate ratio per IQR", "rate difference per IQR"),
                 resid_moran = c(zmi$estimate[1], NA, NA))
print(tz, digits = 4); write.csv(tz, file.path(TAB, "table6_zip_level_asthma.csv"), row.names = FALSE)

hr("10. Pollutant-sign diagnostics")
for (oc in c("asthma", "ami")) for (pv in c("pm25", "ozone")) {
  r <- fit_ols(as.formula(paste(oc, "~", pv)), U, term = pv)
  cat(sprintf("  bivariate %s ~ %s (weighted): %.4f, ZIP-cl p %s\n", oc, pv, r$zip["est"], pfmt(r$zip["p"])))
  log_spec("P1", oc, "primary", "WLS bivariate", r$n, pv, r$zip["est"], r$zip["se"], r$zip["lo"], r$zip["hi"], r$zip["p"])
}
cat("\n  Pearson correlations (primary sample):\n")
print(round(cor(U[, c("dist", "pm25", "ozone", "poverty", "density", "asthma", "ami")]), 3))
X <- U[, c("dist", "pm25", "ozone", "poverty")]
vif <- sapply(names(X), function(v) 1 / (1 - summary(lm(X[[v]] ~ ., data = X[, names(X) != v, drop = FALSE]))$r.squared))
cat("  VIF:", paste(names(vif), round(vif, 3), collapse = " | "), "\n")
cat("  Moran's I of each variable (queen):\n")
for (v in c("asthma", "ami", "dist", "pm25", "ozone", "poverty"))
  cat(sprintf("    %-8s %.3f\n", v, moran.test(U[[v]], QU$lw, zero.policy = TRUE)$estimate[1]))

hr("11. Table 9 and the other models in the specification log")
Wm <- as(QU$lw, "CsparseMatrix"); trMC <- trW(Wm, type = "MC")
for (oc in c("asthma", "ami")) {
  f <- if (oc == "asthma") f_a else f_m
  # spatialreg::lagsarlm has no weights argument, so the lag and Durbin models are unweighted
  lag <- lagsarlm(f, data = U, listw = QU$lw, zero.policy = TRUE)
  im <- summary(impacts(lag, tr = trMC, R = 1000), zstats = TRUE, short = TRUE)
  cat(sprintf("\n  [%s] spatial lag (unweighted): rho %.3f\n", oc, lag$rho)); print(im)
  for (k in c("direct", "indirect", "total")) {
    i <- grep("^dist ", rownames(im$semat)); cap <- c(direct = "Direct", indirect = "Indirect", total = "Total")
    est <- im$res[[k]][i]; se <- im$semat[i, cap[[k]]]
    log_spec("A1", oc, "primary", paste("spatial lag (unweighted),", k, "effect"), nrow(U), "dist per mile", est, se, est - 1.96 * se, est + 1.96 * se, im$pzmat[i, cap[[k]]], "simulation-based p-value")
  }
  dur <- lagsarlm(f, data = U, listw = QU$lw, Durbin = TRUE, zero.policy = TRUE)
  imd <- summary(impacts(dur, tr = trMC, R = 1000), zstats = TRUE, short = TRUE)
  cat(sprintf("  [%s] spatial Durbin (unweighted): rho %.3f\n", oc, dur$rho)); print(imd)
  for (k in c("direct", "indirect", "total")) {
    i <- grep("^dist ", rownames(imd$semat)); cap <- c(direct = "Direct", indirect = "Indirect", total = "Total")
    est <- imd$res[[k]][i]; se <- imd$semat[i, cap[[k]]]
    log_spec("A2", oc, "primary", paste("spatial Durbin (unweighted),", k, "effect"), nrow(U), "dist per mile", est, se, est - 1.96 * se, est + 1.96 * se, imd$pzmat[i, cap[[k]]], "simulation-based p-value")
  }
  for (ix in c("poverty", "pm25")) {
    mi <- lm(as.formula(paste(oc, "~ dist *", ix, "+ pm25 + ozone + poverty")), data = U, weights = pop)
    r <- term_row(mi, paste0("dist:", ix), zip_vc(mi, U), length(unique(U$ZIP)) - 1)
    cat(sprintf("  [%s] dist x %s interaction: %.4f (p %s)\n", oc, ix, r["est"], pfmt(r["p"])))
    log_spec("A3", oc, "primary", paste("WLS interaction dist x", ix), nobs(mi), paste0("dist:", ix), r["est"], r["se"], r["lo"], r["hi"], r["p"])
  }
  me <- fit_ols(update(f, . ~ . + I(density / 1000)), U)
  log_spec("A4", oc, "primary", "WLS + density", me$n, "dist per mile", me$zip["est"], me$zip["se"], me$zip["lo"], me$zip["hi"], me$zip["p"])
  Ui <- U[!is.na(U$uninsured), ]
  mu <- fit_ols(update(f, . ~ . + I(density / 1000) + uninsured), Ui)
  log_spec("A5", oc, "primary (uninsured available)", "WLS + density + % uninsured", mu$n, "dist per mile", mu$zip["est"], mu$zip["se"], mu$zip["lo"], mu$zip["hi"], mu$zip["p"])
  cat(sprintf("  [%s] + density: %.3f (p %s) | + density + uninsured: %.3f (p %s, n %d)\n", oc, me$zip["est"], pfmt(me$zip["p"]), mu$zip["est"], pfmt(mu$zip["p"]), mu$n))
  m3 <- lm(update(f, paste(oc, "~ . + factor(substr(ZIP, 1, 3))")), data = U, weights = pop)
  r <- term_row(m3, "dist", zip_vc(m3, U), length(unique(U$ZIP)) - 1)
  log_spec("A6", oc, "primary", "WLS 3-digit ZIP fixed effects", nobs(m3), "dist per mile", r["est"], r["se"], r["lo"], r["hi"], r["p"])
  cat(sprintf("  [%s] 3-digit-ZIP fixed effects: %.3f (p %s)\n", oc, r["est"], pfmt(r["p"])))
}

hr("12. Appendix figures")
XLAB <- "Mean street-network distance to\nthree nearest parks (miles)"
th <- theme_minimal(base_size = 13) + theme(plot.title = element_text(size = 12.5, face = "bold"),
  panel.border = element_rect(colour = "grey70", fill = NA, linewidth = 0.3), plot.margin = margin(6, 10, 6, 6))
sc <- function(yvar, ylab, title) ggplot(U, aes(x = dist, y = .data[[yvar]])) +
  geom_point(colour = "royalblue", alpha = 0.45, size = 0.8) +
  geom_smooth(method = "lm", formula = y ~ x, mapping = aes(weight = pop), colour = "red", fill = "grey55", linewidth = 0.8) +
  labs(x = XLAB, y = ylab, title = title) + th
f7 <- (sc("asthma", "Asthma ED visits per 10,000", "Asthma") + sc("ami", "AMI ED visits per 10,000", "AMI") +
       sc("pm25", "PM2.5 (µg/m³)", "PM2.5") + sc("ozone", "Ozone (ppm)", "Ozone") +
       sc("poverty", "Below twice the federal poverty level (%)", "Poverty")) +
  plot_layout(ncol = 3) + plot_annotation(tag_levels = "A")
ggsave(file.path(FIG, "Figure7.png"), f7, width = 14.5, height = 9.2, dpi = 300, bg = "white")
A <- main$asthma$A$m
f6a <- ggplot(data.frame(fitted = fitted(A), resid = residuals(A)), aes(fitted, resid)) +
  geom_point(shape = 1, colour = "grey25", size = 1, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE, colour = "red", linewidth = 0.7) +
  labs(x = "Fitted values", y = "Residuals", title = "Residuals versus fitted values, adjusted asthma model") + th
cm <- cor(U[, c("dist", "asthma", "ami", "pm25", "ozone", "poverty")])
dimnames(cm) <- rep(list(c("Park distance", "Asthma", "AMI", "PM2.5", "Ozone", "Poverty")), 2)
f6b <- ggcorrplot(cm, type = "upper", lab = TRUE, lab_size = 4.2, outline.color = "grey60",
                  colors = c("#B2182B", "white", "#2166AC")) +
  labs(title = "Pairwise correlations") + theme(plot.title = element_text(size = 12.5, face = "bold"), legend.title = element_blank())
ggsave(file.path(FIG, "Figure6.png"), (f6a + f6b) + plot_annotation(tag_levels = "A"), width = 13.5, height = 5.4, dpi = 300, bg = "white")
cat("  Figures 6 and 7 written to output/figures\n")

hr("13. Tract analysis dataset with sample-membership flags")
ds <- data.frame(geoid_2010 = d$GEOID, zip = d$ZIP, population = d$pop,
                 network_dist_3parks_miles = d$dist, mean_dist_3parks_miles = d$dist_sl,
                 asthma_ed_per10k = d$asthma, ami_ed_per10k = d$ami, pm25_ug_m3 = d$pm25, ozone_ppm = d$ozone,
                 poverty_pct_below_2x_fpl = d$poverty, ozone_ppm_ces3_2012_2014 = d$ozone_ces3,
                 land_area_sq_mi = d$area_sqmi,
                 in_primary = as.integer(d$primary), within_1p25mi = as.integer(d$dist_sl <= 1.25))
check(sum(ds$in_primary) == 2305 && sum(ds$within_1p25mi) == 2291, "sample flags (2,305 primary; 2,291 within 1.25 miles = originally submitted sample)")
write.csv(ds, file.path(OUT, "tract_analysis_dataset.csv"), row.names = FALSE, na = "")
cat("  written: output/tract_analysis_dataset.csv (", nrow(ds), "tracts; missing values left blank )\n")

hr("14. Specification log")
spec <- do.call(rbind, SPEC); rownames(spec) <- NULL
write.csv(spec, file.path(OUT, "spec_log.csv"), row.names = FALSE)
cat("  models/terms logged:", nrow(spec), "\n")
print(sessionInfo())
cat("\nDone\n")
