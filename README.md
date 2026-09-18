# Data and code: park proximity, air pollution, poverty, and ED visit rates in Los Angeles County

These are the tract-level data and the R code behind the Cureus article "Environmental and Socioeconomic Correlates of Emergency Department Visit Rates for Asthma and Acute Myocardial Infarction in Los Angeles County: A Cross-Sectional Ecological Study". The study is a cross-sectional ecological analysis of Los Angeles County census tracts (2010 boundaries). It relates park proximity (the mean street-network distance from each tract's population-weighted center to the three nearest recreational spaces), fine particulate matter (PM2.5), ozone and poverty to age-adjusted emergency department (ED) visit rates for asthma and acute myocardial infarction (AMI).

## Reproducing the analysis
Use R 4.6.1 or later. The results were produced with sf 1.1-2, spdep 1.4-2, spatialreg 1.4-3, sandwich 3.1-3, lmtest 0.9-40, MASS 7.3-65, mgcv 1.9-4, ggplot2 4.0.3, patchwork 1.3.2 and ggcorrplot 0.3.0; `01_prepare_inputs.R` also needs readxl, and `02_network_distance.R` calls Python 3 (3.14 was used; standard library only) for its shortest-path step. The full session information is at the end of `analysis/output/analysis-log.txt`.

From the repository root:

```bash
Rscript analysis/analysis.R > analysis/output/analysis-log.txt 2>&1
```

`analysis.R` runs on the files included here (`data/` and `public-sources/derived/`), so nothing needs to be downloaded. It first reproduces exactly, from `data/oct17csv.csv`, the models of the manuscript as first submitted to the journal (unweighted, with the missing values described below coded as 0); the code, log and output files label these "originally submitted" or "as originally coded". It then runs the analyses reported in the article. It takes one to two minutes and rewrites its outputs in `analysis/output/`.

To rebuild `public-sources/derived/` from the public sources as well, download the files listed in `public-sources/MANIFEST.csv` (URL, retrieval date, size and SHA-256 for each) into `public-sources/downloads/`, saving each under the name in the `file` column. The three `.geojson` entries are ArcGIS feature services, not files: export all features of each layer as GeoJSON (the checksums describe the exports made on the retrieval date, so a new export may not match them byte for byte). Then run, in this order:

```bash
Rscript analysis/01_prepare_inputs.R > analysis/output/01_prepare_inputs.log 2>&1
Rscript analysis/02_network_distance.R > analysis/output/02_network_distance.log 2>&1
```

- `01_prepare_inputs.R` also takes one to two minutes and overwrites `public-sources/derived/`; a rebuilt `zcta_la_2010.gpkg` differs from the deposited one only by its embedded write timestamp.
- `02_network_distance.R` builds the road network from the TIGER/Line edges file, calls `analysis/network_distance.py` for the shortest paths, and writes `public-sources/derived/tract_network_distance.csv` (about one minute). It writes its intermediate network files to a temporary directory; to keep them, pass a directory as its first argument.
- All three scripts stop if any consistency check fails.
- Every statistic in the article comes from `analysis/01_prepare_inputs.R` (the park-layer checks and the recomputed exposure; see `analysis/output/01_prepare_inputs.log`), `analysis/02_network_distance.R` (the street-network distances and walking times; see `analysis/output/02_network_distance.log`) or `analysis/analysis.R` (everything else).
- Every analysis model that `analysis.R` fits is listed in `analysis/output/spec_log.csv`. Auxiliary fits used only to compute partial R², variance inflation factors and the density-interaction Wald test, and the fitted lines in Figures 6 and 7, are not listed separately. The full printed output is in `analysis/output/analysis-log.txt`.
- Three implementation details needed to reproduce the printed digits exactly. The six-nearest-neighbor weights are built from the tract centroids in `tract_exposures.csv` (`centroid_x`, `centroid_y`). In the ZIP-level analysis, neighbors are queen-contiguous areas among the 259 analyzed, the one area with no contiguous neighbor (90704, Santa Catalina Island) is linked to its nearest area, and Moran's I uses Pearson residuals. The spatial lag and Durbin estimates in Table 9 use `spatialreg::impacts` with Monte Carlo traces (`trW(type = "MC")`) and 1,000 simulations; their CIs are the estimate ± 1.96 simulated standard errors, and their p-values come from the simulated distribution. The random seed is set once, at the top of `analysis.R` (`set.seed(20260911)`), so these simulated values also depend on the random draws made earlier in the script (the scrambled-residual check and the Monte Carlo traces).
- The `poly2nb` and `knn2nb` warnings in the log are expected. The contiguity graph of the primary sample has three disconnected parts and no tract without a neighbor. In some restricted samples a few tracts have no contiguous neighbor, and those models are fitted with `zero.policy = TRUE`. The "65 sub-graphs" message comes from the one-nearest-neighbor step that is used only to find the nearest area for ZIP code 90704, not from the weights used in any model.
- p-values written as 0 or 2.2e-16 in the log and the CSV files are below the numerical precision of the routine that produced them; read them as < 1e-15.
- The maps in the article (Figures 1-5) were drawn separately from the same data and are not produced by these scripts.

## Contents
```
data/                                  the author's original data files, unmodified
  oct17csv.csv                         merged analysis file exported from ArcGIS (2,343 tracts x 73 columns)
  Outlier_Tracts_FINAL.csv             the 52 tracts with mean straight-line park distance > 1.25 miles (tract IDs as 11-character text; values rounded to two decimals, so two rows show 1.25; not read by the scripts)
  revisedLAcensustracts.zip            tract layer used for centroids and joins (2,343 tracts, 2010 boundaries)
  revisedLAcensustractslayer.zip       the same tract layer as republished in ArcGIS Online (Web Mercator, with the joined MEAN_Total field; not read by the scripts)
public-sources/
  MANIFEST.csv                         public source files: URL, retrieval date, size, SHA-256
  derived/                             inputs built by 01_prepare_inputs.R and 02_network_distance.R (used by analysis.R)
  downloads/                           placeholder only; put the MANIFEST files here to run 01_prepare_inputs.R and 02_network_distance.R
analysis/
  01_prepare_inputs.R                  recomputes park-access metrics and builds the tract- and ZIP-level inputs
  02_network_distance.R                street-network distances and walking times from population-weighted tract centers
  network_distance.py                  shortest-path step called by 02_network_distance.R (Python 3, standard library)
  analysis.R                           all statistical analyses, tables and appendix figures
  output/                              logs, tables (CSV), specification log, figures, analysis dataset
GIS-COMMAND-SEQUENCE.md                how the park-distance exposure was built in ArcGIS, and how the street-network distance was measured
LICENSE                                MIT License (code)
```

**Output files and the article.** The numbers in the table file names are the script's own; from Table 5 onward they differ from the article's table numbers, which the log headings use.

| File in `analysis/output/` | In the article |
|---|---|
| `tables/table1_descriptives.csv` | Table 1 |
| `tables/table2_within_vs_beyond_1p25mi.csv` | Table 2 |
| `tables/table3_bivariate.csv` | Table 3 |
| `tables/table4a_main_models.csv`, `tables/table4a_summary.csv` | Table 4 |
| `tables/table4b_within_between.csv`, `tables/table4b_variance_shares.csv`; `spec_log.csv`, rows T4B (ZIP-code fixed effects) | Table 5 (in `table4b_within_between.csv` and the T4B rows of `spec_log.csv`, `lo` and `hi` are per unit of the variable, that is per mile for park distance; Table 5 reports park distance per IQR, the per-mile value multiplied by `iqr_miles` in `table4a_summary.csv`. `table4b_variance_shares.csv` holds only the variance shares and has no confidence limits) |
| `tables/table5_sensitivity.csv` | Tables 6 and 7 (its `pov_*` columns hold each sensitivity model's poverty coefficient, which the article does not tabulate) |
| `tables/table6_zip_level_asthma.csv` | Table 8 |
| `spec_log.csv`, rows A1-A6 | Table 9 |
| `figures/Figure6.png`, `figures/Figure7.png` | Figures 6 and 7 |

## Missing values that appear as 0 in the files in `data/`
- In the official CalEnviroScreen 4.0 release, some tract values are missing. The shapefile download from the Office of Environmental Health Hazard Assessment (OEHHA) codes them -999 and the official spreadsheet shows them as NA; the tract layer in `data/` and `data/oct17csv.csv` carry them as 0. They are:
  - poverty in 38 tracts: the 16 tracts with no residents, plus 22 tracts whose poverty estimate OEHHA judged unreliable;
  - asthma and AMI rates in 9 of the zero-resident tracts;
  - the two age percentages used here (under 10, 65 and over) in the 16 zero-resident tracts. In one further tract (06037930401) the tract layer and the export carry these values as 0: OEHHA's shapefile download codes them -999, although the Demographic Profile sheet of the official spreadsheet reports values for it (10.2 and 21.06). The spreadsheet values were used. All other zeros in these columns are true zeros in the official data.
- Population, PM2.5, ozone, asthma, AMI and poverty in `oct17csv.csv` otherwise match the official spreadsheet exactly, and the two age percentages agree with it to two decimals (the export carries the OEHHA shapefile's four-decimal values, which the spreadsheet rounds). `analysis.R` checks each of these, value by value, against `public-sources/derived/ces4_official_la.csv`. The other columns of `oct17csv.csv` (the remaining CalEnviroScreen 4.0 fields, with names shortened by the shapefile format, and ArcGIS bookkeeping fields) are not used; OEHHA describes the CalEnviroScreen fields in the Data Dictionary sheet of the workbook listed in `public-sources/MANIFEST.csv`.
- The analysis reported in the article treats these values as missing. The primary analysis uses the 2,305 tracts with residents and complete data.
- The originally submitted models, which used the zeros and the straight-line park-distance measure, are reproduced as a check and reported as sensitivity analyses (the three "Straight-line measure ... missing values coded as zero" rows in each of Tables 6 and 7).

## Data dictionary
**Study variables** in `data/oct17csv.csv`. Tract values are indicator values, not percentiles; for example, `Asthma` is used, not `Asthma_Pct`.

| Column | Meaning | Source |
|---|---|---|
| `Tract` | Census tract GEOID stored as a number (10 digits, no leading zero) | CalEnviroScreen 4.0 |
| `ZIP` | ZIP code assigned to the tract by CalEnviroScreen (used for ZIP-clustered standard errors) | CalEnviroScreen 4.0 |
| `Population` | Tract population (American Community Survey, ACS, 2015-2019) | CalEnviroScreen 4.0 |
| `MEAN_Total_Miles` | The straight-line park-distance measure of the initial analysis (miles): mean straight-line distance from the tract centroid to the nearest edge of each of the three nearest recreational spaces. The article reports it as a sensitivity analysis ("straight-line measure"); the primary exposure is `net3_popcentroid` in `public-sources/derived/tract_network_distance.csv` | Derived (see `GIS-COMMAND-SEQUENCE.md`) |
| `Asthma` | Age-adjusted asthma ED visits per 10,000 residents, 2015-2017 (spatially modeled) | CalEnviroScreen 4.0 |
| `Cardiovasc` | Age-adjusted AMI ED visits per 10,000 residents, 2015-2017 (spatially modeled) | CalEnviroScreen 4.0 |
| `PM2_5` | Annual mean PM2.5 (µg/m³), 2015-2017 | CalEnviroScreen 4.0 |
| `Ozone` | Mean of summer months (May-October) of the daily maximum 8-hour ozone concentration (ppm), 2017-2019 | CalEnviroScreen 4.0 |
| `Poverty` | % of residents below twice the federal poverty level (ACS 2015-2019) | CalEnviroScreen 4.0 |
| `Shape_Area` | Tract area in m² (California Albers). Use this, not `Shape__Area`, which is in Web Mercator and inflated | ArcGIS export |
| `Child_10`, `Elderly_65` | % under 10 and % 65 and over (ACS 2015-2019 estimates), used in the ZIP-code count models | CalEnviroScreen 4.0 |
| `FREQUENCY` | Number of park distances averaged (3 on every row) | ArcGIS Summary Statistics |

**How the ED outcomes were produced.** From the CalEnviroScreen 4.0 report, pp. 152-153 and 157-158:
1. ZIP-code ED visit rates were age-adjusted.
2. They were "spatially modeled to provide estimates for ZIP codes with fewer than 12 ED visits."
3. ZIP values were assigned to 2010 census blocks by areal apportionment.
4. Tract rates are population-weighted averages of those block values.

**Derived inputs** (`public-sources/derived/`):

| File | Contents |
|---|---|
| `ces4_official_la.csv` | Official CalEnviroScreen 4.0 values for the 2,343 tracts, missing values written as `NA`: `GEOID`, `population`, `asthma`, `ami`, `poverty`, `pm25`, `ozone`, and from the workbook's Demographic Profile sheet `child_10`, `elderly_65` (percentages) |
| `ces3_ozone_la.csv` | CalEnviroScreen 3.0 ozone, 2012-2014 (same metric as 4.0): `GEOID`, `ozone_ces3` |
| `tract_exposures.csv` | Park-access metrics recomputed from geometry (columns below) |
| `tract_network_distance.csv` | Street-network distances and walking times from the population-weighted tract centers, built by `02_network_distance.R` (columns below) |
| `zcta_tract_rel_la.csv` | 2010 ZIP Code Tabulation Area (ZCTA)-to-tract relationship records for Los Angeles County (Census Bureau): `ZCTA5`, `GEOID`, `POPPT` (2010 population in the overlap), `ZPOP` (2010 ZCTA population), `ZPOPPCT` (share of the ZCTA population in the overlap, %) |
| `cdph_asthma_zip_la.csv` | Asthma ED visits by ZIP code for 2015-2017, aggregated for this study from the California Department of Public Health table (all ages): `zip`; `years` (number of years with published values); `visits_2015_2017` (sum of the annual counts); `mean_age_adj_rate` (mean of the annual age-adjusted rates per 10,000 residents). Counts below 12 are not published |
| `zcta_la_2010.gpkg` | 2010 ZCTA polygons for Los Angeles County (Census TIGER/Line), California Albers (EPSG:3310) |
| `acs_uninsured_2015_2019.csv` | `Tract`, `pct_uninsured`: % uninsured (ACS 2015-2019), used in one extended model. Built by `01_prepare_inputs.R` from the Census Bureau's ACS 2015-2019 5-year summary file for California (sequence file `20195ca0131000.zip`, table B27001, with the geography file `g20195ca.csv`; both in `MANIFEST.csv`): the sum of the 18 "No health insurance coverage" cells divided by the total (B27001_001) × 100. It has the 2,323 Los Angeles County tracts with a nonzero denominator; 21 study tracts have no value (one is not in the 2019 tract list and 20 have a zero denominator), one of them in the primary sample, so the model using it has n = 2,304 |

**Columns of `tract_exposures.csv`** (all distances in miles):

| Column | Meaning |
|---|---|
| `GEOID` | Tract GEOID (11 characters) |
| `recorded_dist3_miles` | The straight-line measure (= `MEAN_Total_Miles`) |
| `centroid_x`, `centroid_y` | Tract centroid placed as ArcGIS Find Centroids does, California Albers (EPSG:3310), meters |
| `popcentroid_x`, `popcentroid_y` | Census 2010 population-weighted tract center, EPSG:3310, meters |
| `dist3_planar`, `dist3_geodesic` | The exposure recomputed from geometry (planar California Albers / spherical geodesic). The "recomputed from geometry" sensitivity row uses whichever of the two correlates better with the recorded values in the primary sample, which is the geodesic version (see the analysis log) |
| `dist1_planar` | Distance to the single nearest park edge |
| `dist3_popcentroid`, `dist1_popcentroid` | Same distances from the Census 2010 population-weighted tract center |
| `acres_0p5mi`, `acres_1mi` | Park acres within 0.5 and 1 mile of the tract centroid (overlaps dissolved) |
| `dist3_pna2016` | Mean distance to the three nearest parks in the 2016 county Park Needs Assessment inventory (excluding "No Public Access") |
| `walkshed2016_share` | Share of the tract's area inside a 2016 half-mile pedestrian-network walkshed of a park |
| `popcentroid_in_walkshed2016` | 1 if the population-weighted center lies inside a 2016 walkshed |

**Columns of `tract_network_distance.csv`** (distances in miles along the road network; see `GIS-COMMAND-SEQUENCE.md` for the method):

| Column | Meaning |
|---|---|
| `GEOID` | Tract GEOID (11 characters) |
| `net3_popcentroid` | The article's primary exposure: mean street-network distance from the Census 2010 population-weighted tract center to the three nearest parks |
| `net1_popcentroid` | Street-network distance from that center to the single nearest park |
| `walk3_popcentroid_min`, `walk1_popcentroid_min` | The same two distances as walking time in minutes at 3 miles per hour (20 minutes per mile) |
| `net3_popcentroid_open` | The three-park network distance with the park set restricted to the 2,230 open-access features (sensitivity analysis) |
| `net3_popcentroid_merged` | The three-park network distance with overlapping park features merged into one park each (2,115 parks; sensitivity analysis) |
| `net3_popcentroid_pna2016` | The three-park network distance to the 2016 county Park Needs Assessment inventory (2,828 features, excluding "No Public Access"), so that parks, roads and outcomes fall in one period (sensitivity analysis) |
| `net_snap_m` | Straight-line distance in meters from the population-weighted center to the nearest point of the road network (included in the distances) |
| `net_component` | Part of the road network used: 1 = the main network, 2 = Santa Catalina Island |
| `net3_parks` | Row numbers (1 to 2,509) of the three nearest parks in `parks_recreational_spaces.geojson` (a source download; see `public-sources/MANIFEST.csv`), nearest first |

**Analysis dataset** (`analysis/output/tract_analysis_dataset.csv`, one row per tract, missing values left blank):

| Column | Meaning |
|---|---|
| `geoid_2010`, `zip` | Tract GEOID (11 characters) and its CalEnviroScreen ZIP code |
| `population` | Tract population |
| `network_dist_3parks_miles` | The primary exposure: mean street-network distance from the population-weighted tract center to the three nearest parks |
| `mean_dist_3parks_miles` | The straight-line measure from the tract centroid (the exposure of the initial analysis) |
| `asthma_ed_per10k`, `ami_ed_per10k` | Outcomes, official values (blank where CalEnviroScreen reports none) |
| `pm25_ug_m3`, `ozone_ppm`, `poverty_pct_below_2x_fpl` | Covariates, official values |
| `ozone_ppm_ces3_2012_2014` | CalEnviroScreen 3.0 ozone |
| `land_area_sq_mi` | Tract polygon area in square miles (from `Shape_Area`) |
| `in_primary` | 1 = in the primary sample (residents and complete data; n = 2,305) |
| `within_1p25mi` | 1 = mean straight-line distance ≤ 1.25 miles. This is the originally submitted analytic sample (n = 2,291) |

**Specification log** (`analysis/output/spec_log.csv`): one row per reported model term, with `id` (the part of the analysis), `outcome`, `sample`, `model`, `n`, `term`, `est`, `se`, `lo` and `hi` (95% CI), `p` and `note`. For most rows the unit is given in `term` (for example "dist per mile", "original straight-line dist per mile" or "exposure per 0.376"; `dist` is the street-network exposure unless the term says otherwise). Where `term` is a variable name, the estimate is per unit of that variable (ozone is per 0.01 ppm where `note` says so, that is in the T4 rows, and per 1 ppm in the P1 rows); T4std rows are standardized coefficients, and S7int rows give only the joint-test p-value. The `id` values are: G1, the originally submitted models reproduced; C9, the same models within 1.0 mile; C8, the distance-poverty association under three sample definitions; T3, Table 3; T3u, the same bivariate models unweighted (also in `tables/table3_bivariate.csv`; not shown in Table 3); T3sl, the same weighted bivariate models with the straight-line measure (quoted in the Results; not tabulated); T4 and T4std, Table 4 (the unweighted T4 rows are the per-mile form of the "Primary sample, unweighted" rows of Tables 6 and 7); T4B, Table 5; the S ids, the rows of Tables 6 and 7, in the order of the table (the ids present are S0 to S6, S7T1 to S7T3 for the three density tertiles, S8 to S19 and S22 to S27; there is no plain S7, S20 or S21) (S0 is the primary street-network exposure; S24 the straight-line measure; S1, S2 and S23 the originally submitted samples with the straight-line measure; S22 the street-network distance to the single nearest park; S25 and S26 the open-access-only and merged-polygon variants; S27 the street-network distance to the 2016 county park inventory; S7int is the density-interaction test quoted in the Results); Z1 to Z3, Table 8; P1, bivariate models of each outcome on each pollutant (section 10 of the log); A1 to A6, Table 9. In the table files, columns beginning `ols_` refer to the least-squares model and columns beginning `sem_` to the spatial error model, `lo` and `hi` are 95% confidence limits, and `scale` in `table5_sensitivity.csv` is the increment of the exposure that the coefficients refer to (its IQR in the primary sample, which is also the increment used for the rows fitted on other samples; 1 for log2 distance).

## Sources
- CalEnviroScreen 4.0 and 3.0: California Office of Environmental Health Hazard Assessment (Sacramento, CA, USA).
- Recreational Spaces: Los Angeles County Department of Public Health (Los Angeles, CA, USA).
- 2016 Park Needs Assessment inventory and walksheds: County of Los Angeles Department of Parks and Recreation.
- ZIP-code asthma ED visits: California Department of Public Health, "Asthma ED Visit Rates by ZIP Code 2013-Present", CHHS Open Data Portal (https://data.chhs.ca.gov/dataset/asthma-emergency-department-visit-rates; resource last updated 5 December 2025). `public-sources/derived/cdph_asthma_zip_la.csv` was aggregated from that table for this study; it has been modified from the original and is not official government data.
- Centers of population, ZCTA relationships, ZCTA polygons, ACS estimates and the TIGER/Line 2016 road network (All Lines edges, Los Angeles County): U.S. Census Bureau (Washington, DC, USA). The ZCTA files in `public-sources/derived/` are subsets of the Census Bureau files, reprojected and repackaged for this study; `tract_network_distance.csv` was computed from the road network for this study.

URLs and checksums for the downloaded sources are in `public-sources/MANIFEST.csv`.

## Reuse
The code (R and Python) is released under the MIT License (see `LICENSE`). The values created for this study (the park-distance and park-access measures, the analysis dataset and the output tables) are released under the Creative Commons Attribution 4.0 International License (https://creativecommons.org/licenses/by/4.0/); please cite the article when reusing them. The source data come from the agencies listed above and remain subject to their terms; please cite those agencies when reusing them.

## Privacy
Every file is a publicly available, area-level aggregate. The smallest unit is a census tract or ZIP code. There are no individual records and no personal identifiers.
