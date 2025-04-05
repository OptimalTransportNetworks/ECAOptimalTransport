library(fastverse)
fastverse_extend(sf, haven)

# Load DTA
ECA_NUTS <- read_dta("data/ECA_database_country and nuts_apc_unique.dta") # |> funique()
# ECA_NUTS |> write_dta("data/ECA_database_country and nuts_apc_unique.dta")
fndistinct(ECA_NUTS)
ECA_NUTS |> fselect(year, nuts2id, nuts3id, nuts2idd, nuts3idd) |> fnunique()
# ECA_NUTS |> gvr("year|nuts") |> fnunique()

# Load geometries
ECA_shp <- haven::read_dta("data/ECA_shp/ECA_shp.dta") |> rm_stub("_")
fndistinct(ECA_shp)

# Load GADM
st_layers("/Users/sebastiankrantz/Documents/Data/GADM/gadm_410-levels.gpkg")
GADM2 <- st_read("/Users/sebastiankrantz/Documents/Data/GADM/gadm_410-levels.gpkg", layer = "ADM_1") 

# Generate bbox from ECA shp
ECA_shp_sf <- ECA_shp |> 
  fsubset(is.finite(X) & is.finite(Y) & between(X, -180, 180) & between(Y, -90, 90)) |> 
  st_as_sf(coords = c("X", "Y"), crs = 4326) 

ECA_bbox <- ECA_shp_sf |> 
  st_bbox() |> 
  st_as_sfc() |>
  st_buffer(dist = units::as_units(500, "km")) 

st_crs(ECA_bbox) <- st_crs(GADM2)

# Subset GADM2 to ECA bbox
qtable(st_is_valid(GADM2$geom))
GADM2 %<>% subset(st_is_valid(geom)) %>% subset(st_intersects(geom, ECA_bbox, sparse = FALSE))
# mapview::mapview(st_as_sf(ECA_bbox)) + mapview::mapview(GADM2)

# Map firms to administrative areas
map <- st_intersects(ECA_shp_sf, GADM2)
qtable(vlengths(map))
# ECA_shp_sf |> subset(vlengths(map) == 0) |> mapview::mapview()
for (i in c(1, 5, 10)) { #
  map2 <- st_within(ECA_shp_sf, st_buffer(GADM2, dist = units::as_units(i, "km")))
  # qtable(vlengths(map2)[vlengths(map) == 0])
  map[vlengths(map) == 0] <- map2[vlengths(map) == 0]
}; rm(map2)
qtable(vlengths(map))
ECA_shp_sf$GADM2 <- ffirst(map) # fsum(ffirst(map, drop = FALSE))
GADM2 %<>% ss(funique(ffirst(map), sort = TRUE))

# Loading NUTS3
NUTS3 <- readRDS("data/NUTS3_pop.rds") |> 
  st_transform(4326) |>
  subset(levl_code == 3 & 
         cntr_code %!in% c("UK", "IE", "IS", "NO", "SE", "FI", "CY") & 
         Reduce("&", Map(`>`, mctl(st_coordinates(st_centroid(geometry))), list(-10, 28))))

# NUTS3 <- readRDS("data/NUTS3_pop.rds")
setdiff(ECA_NUTS$nuts3id, NUTS3$nuts_id)
setdiff(NUTS3$nuts_id, ECA_NUTS$nuts3id)

# Pretty good overlap...
# mapview::mapview(NUTS3) + mapview::mapview(GADM2)

# Combining 
# mapview::mapview(NUTS3) + mapview::mapview(subset(GADM2, GID_0 %!in% c("POL", "HUN", "SVN", "SVK", "ALB", "ROU", "HRV", "SRB", "BGR", "MNE", "MKD", "GRC", "TUR")))

NUTS3_ext <- NUTS3 |> st_cast("MULTIPOLYGON") |>
  rowbind(GADM2 |> 
      subset(GID_0 %!in% c("POL", "HUN", "SVN", "SVK", "ALB", "ROU", "HRV", "SRB", "BGR", "MNE", "MKD", "GRC", "TUR")) |>
      st_make_valid() |>
      rmapshaper::ms_simplify(keep = 0.3) |>
      st_cast("MULTIPOLYGON") |>
      unclass() |>
      fcompute(leve_code = 3L,
               nuts_id = GID_1, 
               cntr_code = iif(GID_0 == "XKO", "KO", countrycode::countrycode(GID_0, "iso3c", "iso2c")),
               name_latn = VARNAME_1, 
               nuts_name = NAME_1, 
               geometry = geom), 
  fill = TRUE) 

NUTS3_ext %<>% st_make_valid() %>% 
  rmapshaper::ms_simplify(keep = 0.5) %>%
  st_make_valid() %>% 
  colorderv(is_categorical)
  
# mapview::mapview(NUTS3_ext[, c("nuts_id", "nuts_name")])

# Interior Polygons. Example: Ansbach, Germany
NUTS3_ext %<>% sfheaders::sf_remove_holes() 

st_within(subset(NUTS3_ext, nuts_id == "DE251", geometry),
          subset(NUTS3_ext, nuts_id == "DE256", geometry)) 

# Now removing interior polygons
within_matrix <- st_within(NUTS3_ext, NUTS3_ext, sparse = FALSE)
diag(within_matrix) <- FALSE
interior_ids <- which(rowSums(within_matrix) > 0)
# mapview::mapview(NUTS3_ext[interior_ids, c("nuts_id", "nuts_name")])
container_ids <- apply(within_matrix, 1, function(x) which(x)[1])
container_ids <- container_ids[interior_ids]

NUTS3_ext[interior_ids, cat_vars(NUTS3_ext, "names")] <- 
  NUTS3_ext[container_ids, cat_vars(NUTS3_ext, "names")]

NUTS3_ext %<>% collap(~ nuts_id, fsum, ffirst)

# Getting overall Boundary and Nightlights
NUTS3_ext_hull <- st_union(NUTS3_ext) |> st_convex_hull() |> st_as_sf()
# mapview::mapview(NUTS3_ext_hull)

# bearer <- "eyJ0eXAiOiJKV1QiLCJvcmlnaW4iOiJFYXJ0aGRhdGEgTG9naW4iLCJzaWciOiJlZGxqd3RwdWJrZXlfb3BzIiwiYWxnIjoiUlMyNTYifQ.eyJ0eXBlIjoiVXNlciIsInVpZCI6InNlYmtyYW50eiIsImV4cCI6MTc0ODMwNzE5NCwiaWF0IjoxNzQzMTIzMTk0LCJpc3MiOiJodHRwczovL3Vycy5lYXJ0aGRhdGEubmFzYS5nb3YiLCJpZGVudGl0eV9wcm92aWRlciI6ImVkbF9vcHMiLCJhY3IiOiJlZGwiLCJhc3N1cmFuY2VfbGV2ZWwiOjN9.CTnC0Do3llf-wt_hba-OkGJBhXzKRshWBFpK4cx3xtJaLauLS7SiAeUrMR1WqAbq9h4R30E8PUCVnnHRMZd5uf5D6VTXxhhqEVT-4yUt33Ji1Q6pzHUwIb30uqKMkggxb_n4gL-8s3oot4b6PwutiU4bgsZX_YALakTeLHnLg8Ry9r9zVFmOJeorYgCpWRsx0XdqT-_YNmYA8MrK9RNZ3FmFh2Hs3_xN_KPBgMGEb6jpNciBQiir1cjvFQuBPWGP-Eix1PfmBN-DJuKXdp35tKGduJtUTpKxJ-n_ZJLCgJZ5Jw63UewHJduu0kziPy3IYBO8ufMkufv3vv_uSyvK7A"
# NL23 <- blackmarbler::bm_raster(roi_sf = NUTS3_ext_hull,
#                   product_id = "VNP46A4",
#                   date = "2023",
#                   bearer = bearer) |>
#   terra::aggregate(fact = 4, fun = "sum")
# NL23 <- NL23 |> as.data.frame(xy = TRUE) |> fsubset(t2023 > 0)
# qs::qsave(NL23, "data/NL23_ECA_agg.qs")
NL23 <- qs::qread("data/NL23_ECA_agg.qs")

# Find activity centroids
m <- st_within(st_as_sf(NL23, coords = c("x", "y"), crs = 4326), NUTS3_ext)
NL23 %<>% fsubset(vlengths(m) > 0)
NL23$nuts_row <- as.integer(ffirst(m[vlengths(m) > 0]))
rm(m)
NL23 %<>% collap( ~ nuts_row, w = ~ t2023)
# Inspect
qsu(NL23)
# NL23 |> st_as_sf(coords = c("x", "y"), crs = 4326) |> 
#   mapview::mapview(zcol = "t2023", layer.name = "Nightlights 2023") + 
#   mapview::mapview(NUTS3_ext)
# # Removing shapes for which there are no (non-zero) nightlights (Estonia north)
# NUTS3_ext |> ss(setdiff(seq_row(NUTS3_ext), NL23$nuts_row)) |> mapview::mapview()
NUTS3_ext %<>% ss(NL23$nuts_row)

# Centroids not in the shape
not_inside <- which(!s2::s2_within(with(NL23, s2::s2_lnglat(x, y)), NUTS3_ext$geometry))

# Getting nearest point
nearest_points <- st_nearest_points(NUTS3_ext[not_inside, ],
  NL23[not_inside, ] |> st_as_sf(coords = c("x", "y"), crs = 4326), pairwise = TRUE) |> 
  st_cast("POINT", group_or_split = FALSE, warn = FALSE) |> st_as_sf()
  
# NL23[not_inside, ] |> 
#   st_as_sf(coords = c("x", "y"), crs = 4326) |> 
#   mapview::mapview(zcol = "t2023", layer.name = "Nightlights 2023") + 
#   mapview::mapview(NUTS3_ext[not_inside, ]) +
#   mapview::mapview(nearest_points)

# Replacing 
set(NL23, not_inside, c("x", "y"), mctl(st_coordinates(nearest_points)))

tfm(NUTS3_ext) <- frename(NL23, x = lon, y = lat, t2023 = nightlights)
saveRDS(NUTS3_ext, "data/NUTS3_ext.rds")

# Now Computing Distance Matrix
source("code/helpers.R")
library(osrm)

dist <- split_large_dist_matrix(qDF(qM(fselect(qDF(NUTS3_ext), nuts_id, lon, lat), 1)))
anyNA(dist$distances)
anyNA(dist$durations)

# Road Transport Cost Estimation
# ChatGPT suggests: Cost [€] = 1.05 × (distance in km) + 0.60 × (travel time in minutes)
# Iimi (2023) suggests: log(U.S. cents per ton-km) = 4.650 - 0.395 × log(speed in km/h) - 0.064 × log(distance in km) + 0.024 × crossborder dummy
# Also sector-specific results and different equations for domestic/international shipments...

# Compute travel cost
dist$travel_cost <- dist$durations / 3600









# # Fetch nightlights for NUTS1 regions
# NUTS1 <- NUTS |> subset(levl_code == 1) |> st_transform(4326)
# #### Define NASA bearer token
# bearer <- "eyJ0eXAiOiJKV1QiLCJvcmlnaW4iOiJFYXJ0aGRhdGEgTG9naW4iLCJzaWciOiJlZGxqd3RwdWJrZXlfb3BzIiwiYWxnIjoiUlMyNTYifQ.eyJ0eXBlIjoiVXNlciIsInVpZCI6InNlYmtyYW50eiIsImV4cCI6MTc0ODMwNzE5NCwiaWF0IjoxNzQzMTIzMTk0LCJpc3MiOiJodHRwczovL3Vycy5lYXJ0aGRhdGEubmFzYS5nb3YiLCJpZGVudGl0eV9wcm92aWRlciI6ImVkbF9vcHMiLCJhY3IiOiJlZGwiLCJhc3N1cmFuY2VfbGV2ZWwiOjN9.CTnC0Do3llf-wt_hba-OkGJBhXzKRshWBFpK4cx3xtJaLauLS7SiAeUrMR1WqAbq9h4R30E8PUCVnnHRMZd5uf5D6VTXxhhqEVT-4yUt33Ji1Q6pzHUwIb30uqKMkggxb_n4gL-8s3oot4b6PwutiU4bgsZX_YALakTeLHnLg8Ry9r9zVFmOJeorYgCpWRsx0XdqT-_YNmYA8MrK9RNZ3FmFh2Hs3_xN_KPBgMGEb6jpNciBQiir1cjvFQuBPWGP-Eix1PfmBN-DJuKXdp35tKGduJtUTpKxJ-n_ZJLCgJZ5Jw63UewHJduu0kziPy3IYBO8ufMkufv3vv_uSyvK7A"
# NL21 <- vector("list", length = nrow(NUTS1))
# for (i in seq_along(NL21)[-(1:124)]) {
#   print(NUTS1$nuts_name[i])
#   NL21[[i]] <- bm_raster(roi_sf = NUTS1[i,1],
#                          product_id = "VNP46A4",
#                          date = "2021",
#                          bearer = bearer) |>
#                as.data.frame(xy = TRUE)
# }
# names(NL21) <- NUTS1$nuts_id
# qs::qsave(NL21, "data/NL21.qs")
# 
# NL21 <- rowbind(NL21, idcol = "NUTS1")
