# GIS command sequence: how the park-distance measures were built

The article's primary exposure is the mean street-network distance from each tract's population-weighted center to its three nearest recreational spaces; how it is computed is described in the last section of this file. The first part of the file documents how `MEAN_Total_Miles` in `data/oct17csv.csv` was produced: the mean straight-line distance from each census tract centroid to its three nearest recreational spaces, which the article reports as a sensitivity analysis ("straight-line measure") and uses for Figures 5B and 5C and the 1.25-mile comparisons. Statements are labeled by source:

- [D] documented in the author's own workflow notes, written when the work was done;
- [V] verified empirically from the geometry by `analysis/01_prepare_inputs.R` (output in `analysis/output/01_prepare_inputs.log`).

Software: ArcGIS Online and ArcGIS Pro (Esri, Redlands, CA, USA). The author's file dates place the GIS work in fall 2025 (CalEnviroScreen data folder 16 September 2025, ArcGIS project data 2 October 2025, exported analysis file `oct17csv.csv` 17 October 2025); these dates cannot be checked from the deposited files.

## Inputs
| Layer | Source | Notes |
|---|---|---|
| Census tracts | CalEnviroScreen 4.0 geodatabase/shapefile (OEHHA, Sacramento, CA, USA), https://oehha.ca.gov/calenviroscreen/maps-data | 2010 tract boundaries. Los Angeles County selected with `tract >= 6037000000 AND tract < 6038000000` (2,343 tracts) and exported as a new layer, published to ArcGIS Online as `revisedLAcensustracts` [D] (the same layer is packaged here as `data/revisedLAcensustracts.zip`) |
| Recreational spaces | Los Angeles County Department of Public Health, "Recreational Spaces", Enterprise GIS Hub item `4e6083c7960d4a1caedb57aa76473e4f`; service https://services.arcgis.com/RmCCgQtiZLDCtblq/arcgis/rest/services/Recreational_Spaces/FeatureServer/0 | 2,509 polygons; portal item created 13 January 2024; the service reports its data were last edited 13 January 2024. All features were used, whatever their access type (Open 2,230; Restricted 236; Unknown 43) [V] |

## Steps
1. **Find Centroids** (ArcGIS Online), with input `revisedLAcensustracts` and the option "Contained within each input feature". The output is `CentroidsLA`, one point per tract [D]. Where the geometric centroid falls outside its tract, ArcGIS places the point inside the tract; this affects 16 tracts [V].
2. **Find Closest** (ArcGIS Online) [D]:
   - input layer `CentroidsLA`; near layer Recreational Spaces;
   - measurement type: line (straight-line) distance;
   - number of closest locations: 3 per input;
   - search range limited to 30 miles; units: miles.

   The output is `Nearest Parks LA`, with connecting lines and their straight-line distances. Distances run to the nearest edge of each park polygon, not to the park's centroid [V]: distances to park edges reproduce the recorded values (r = 0.996), whereas distances to park centroids do not (r = 0.891, with 0.6% of tracts within 0.01 miles).
3. **Summary Statistics** (ArcGIS Pro) [D]:
   - input table: the connecting lines;
   - statistics field: straight-line distance (miles); statistic: Mean;
   - case field: the centroid's tract FIPS code.

   This gives one mean per tract (`MEAN_Total_Miles`). `FREQUENCY = 3` on every row confirms that three parks were averaged.
4. **Join** (ArcGIS Pro) [D]:
   - a text tract ID was created with the Field Calculator expression `str(int(!tract!)).zfill(11)`;
   - Add Join linked the summary table to the CalEnviroScreen tract layer;
   - the result was exported as a shapefile and then as `oct17csv.csv`.

## Rules and verification

### Overlapping or zero-area polygons
The layer has no zero-area and no duplicate polygons. 473 of the 2,509 features partly overlap at least one other feature (478 when features lying wholly inside another feature, or wholly containing one, are also counted), because the layer compiles several source datasets. Find Closest treats every polygon feature as a separate location, so overlapping features can each count among a tract's three nearest; no dissolve was applied, and the recomputation below, which also does not dissolve, reproduces the recorded values [V]. The park-area sensitivity measure (acres within 0.5 and 1 mile) dissolves overlaps first, so no area is counted twice.

### Reproduction
`01_prepare_inputs.R` recomputes the three-park mean edge distance from the geometry, using the same centroid rule and the copy of the same service listed in `public-sources/MANIFEST.csv` (data last edited January 2024). It agrees with the recorded values at r = 0.996 (planar version; the geodesic version is almost identical, see below; median absolute difference 0.0011 miles; 89.9% of tracts within 0.01 miles and 91.5% within 0.02 miles). The remaining differences were not traced to a single cause; recomputed from the deposited `public-sources/derived/tract_exposures.csv` (`dist3_planar` against `recorded_dist3_miles`, all 2,343 tracts), 200 tracts differ by more than 0.02 miles, 53 by more than 0.1 mile and 10 by more than 0.25 miles, and the largest difference is 0.76 miles. Models refitted with the recomputed measure give nearly the same estimates as the recorded measure (rows "Straight-line measure recomputed from geometry" in Tables 6 and 7 of the article).

### Geodesic vs planar
The notes record "line distance" but not whether ArcGIS measured it geodesically or on a planar projection. Both recomputations fit the recorded values equally well: planar (California Albers) r = 0.9961 with 89.9% of tracts within 0.01 miles, and geodesic (spherical) r = 0.9961 with 89.8%, so the data cannot distinguish the two.

### Tract vintage
The first attempt used Esri's 2020 "USA Census Tract Boundaries" layer. The work was redone on the 2010 CalEnviroScreen tract layer after the boundaries were found not to match [D]. The recorded distances match centroids of the 2010 layer, as described above, and all 2,343 tract IDs join one-to-one [V].

## Street-network distance (the primary exposure)
`analysis/02_network_distance.R`, with `analysis/network_distance.py` for the shortest paths, measures distance along the road network from each tract's Census 2010 population-weighted center (`cenpop2010_tr06.txt`) to the same 2,509 recreational-space polygons. This is the article's primary exposure (`net3_popcentroid` in `public-sources/derived/tract_network_distance.csv`); the straight-line measure documented above is reported as a sensitivity analysis (the rows labeled "Straight-line measure" in Tables 6 and 7). The log is `analysis/output/02_network_distance.log`.

1. **Network.** TIGER/Line 2016 All Lines edges for Los Angeles County (`tl_2016_06037_edges.zip` in `MANIFEST.csv`; 483,044 edges). Road edges (`ROADFLG = Y`, 399,891) other than primary roads (limited-access highways, MTFCC S1100; 6,122 edges) and ramps (S1630; 9,037 edges) are kept: 384,732 edges, 295,437 nodes (the TIGER node identifiers `TNIDF`/`TNIDT`), 31,291 miles, in California Albers (EPSG:3310) with edge lengths from `sf::st_length`. The network has 120 connected parts; only the parts with at least 500 nodes are used, that is the main network (293,898 nodes, 99.5%) and Santa Catalina Island (744 nodes). The other 118 fragments (at most 341 nodes each) are ignored.
2. **Park access points.** Every used node inside a park polygon or within 60 m of its boundary (34,383 node-park pairs; 123 parks have none), every point where a used road edge crosses a park boundary (7,793 points on 968 parks; no gap), and, for every park, the nearest point on the nearest used road edge. The straight-line distance from the access point to the polygon (0 inside it) is added to the path. The same is done for two sensitivity park sets: the 2,230 open-access features only, and the 2,509 features with overlapping polygons merged (features whose interiors overlap are grouped and unioned; 478 features form 84 merged parks, leaving 2,115 parks); the results are the columns `net3_popcentroid_open` and `net3_popcentroid_merged`. A fourth park set is the county's 2016 Park Needs Assessment inventory (the 2,828 features other than "No Public Access", as used for the straight-line 2016 measure); with it, parks, roads (2016) and outcomes (2015-2017) fall in one period, and the result is the column `net3_popcentroid_pna2016` (r = 0.984 with the primary network measure across all 2,343 tracts; r = 0.970 within the primary sample, which is the value the article reports).
3. **Tract origins.** Each population-weighted center is joined to the nearest point on the nearest used road edge (median 21 m, 95th percentile 85 m; 13 tracts more than 0.25 mile, the largest 2.8 km for mountain tract 06037930301, whose second and third parks are about 28 miles away by road, so its three-park mean is 19 miles) and that straight-line distance is added to the path. A center that lies inside a park or within 60 m of its boundary (170 centers; 74 center-park pairs in which the center lies inside the polygon, which can exceed the number of such centers because park polygons overlap) reaches that park directly at the straight-line distance, the same allowance a road node gets. One center (tract 06037408006) lies 4.5 m from a ten-node fragment that the 500-node rule ignores and is therefore joined to the main network 215 m away; the effect on its distance is at most 0.13 mile.
4. **Shortest paths.** The access points and origins are inserted as nodes into their road edges. One multi-source run of Dijkstra's algorithm, started from all park access points with the park gap as the starting distance, gives every node the distances to its three nearest distinct parks; the mean of the three and the smallest are the tract's `net3_popcentroid` and `net1_popcentroid`. Walking time is the distance at 3 miles per hour.
5. **Checks** (in the log): every used part of the network reaches at least three parks; no network distance is shorter than the straight-line distance from the same center to the same parks; restricting or merging the park set never shortens the three-park distance; the log prints the ratio of network to straight-line distance, the correlations with the straight-line measures, and the number of centers far from the network.

Limitations of the measure, stated in the article: road centerlines are not a pedestrian network (paths, trails and footbridges are missing, while private roads and unpaved vehicular trails are included); every road other than limited-access highways and ramps is treated as walkable; the connectors to the road and to the park boundary are straight lines; the 2016 roads are matched to a 2024 park layer and 2010 population centers.
