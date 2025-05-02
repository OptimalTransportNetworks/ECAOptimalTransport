library(fastverse)
set_collapse(mask = c("manip", "helper", "special"), nthreads = 4)
fastverse_extend(qs, sf, units, stplanr, osrm, install = TRUE)
source("code/helpers.R")
fastverse_conflicts()

# Important Local GDP Data ----------------------------
countries <- wbstats::wb_cachelist$countries %$% iso3c[region_iso3c %==% "ECS"] %>% 
  c("XKO") %>% setdiff(c("GRL", "ISL", "GBR", "IRL", "NOR", "SWE", "FIN", "FRO", "IMN", "CYP"))

# From geojson.io
ECA_outline <- st_read("/Users/sebastiankrantz/Documents/ECAOptimalTransport/data/outline/layers") |> 
               st_buffer(as_units(5, "km"))

# Local GDP data: 0.5 degree, censored for density below 0.01
LGDP <- fread("/Users/sebastiankrantz/Documents/Data/LocalGDP/0_5deg/final_GDP_0_5deg_postadjust_pop_dens_0_01_adjust.csv") |> 
        subset(iso %in% countries) |> 
        slice(cell_id, iso, subcell_id, how = "max", order.by = year) |> 
        transform(longitude = longitude + 0.25, # Coordinate is bottom left corner
                  latitude = latitude + 0.25) %>% 
        collap(~ cell_id + subcell_id, w = ~ pop_cell+1e-5, 
               custom = list(fsum_uw = gvr(., "_GCP_", return = "names"), 
                             fmean = gvr(., "_GDPC_|itude|censored|year|national", return = "names"), 
                             fmode = is_categorical))
  
LGDP %<>% subset(LGDP |> st_as_sf(coords = c("longitude", "latitude"), crs = 4326) |> 
                 st_intersects(ECA_outline, sparse = FALSE))
# shape <- st_read("/Users/sebastiankrantz/Documents/Data/LocalGDP/0_5deg/shapefile") |>
#          base::subset(cell_id %in% LGDP$cell_id)
# shape %<>% join(LGDP, on = c("cell_id", "iso", "subcell_id"), how = "inner", drop = TRUE)
# shape %<>% st_make_valid() 

# Constructing clean shapes
shape <- Map(
  f = function(x, y) {
    st_polygon(list(cbind(x = c(x-0.25, x+0.25, x+0.25, x-0.25, x-0.25), 
                          y = c(y-0.25, y-0.25, y+0.25, y+0.25, y-0.25))))
  },
  x = LGDP$longitude,
  y = LGDP$latitude
) |> st_as_sfc(crs = 4326) |> 
  st_sf(LGDP, geometry = _, crs = 4326)

if(FALSE) mapview::mapview(shape, zcol = "pop_cell") + 
  mapview::mapview(LGDP |> st_as_sf(coords = c("longitude", "latitude"), crs = 4326), zcol = "pop_cell") +
  mapview::mapview(ECA_outline, color = "red")

# Creating Adjacency grid  -----------------------------------------------

# queen = TRUE  →  edge or corner touching
# gap_fix_snap(shape, tol = deg_m(5000)) |> mapview::mapview()
# neigh <- shape |> spdep::poly2nb(queen = TRUE, snap = deg_m(1000))
# table(spdep::card(neigh))
# neigh <- data.frame(from = rep(seq_along(neigh), spdep::card(neigh)),
#                     to   = unlist(unclass(neigh)[spdep::card(neigh) > 0], use.names = FALSE)) 

# Manual: better!
settfmv(shape, c("longitude", "latitude"), round, 2)
neigh <- shape |> with(cbind(longitude, latitude)) |> mrtl() |> 
  lapply(function(r) {
    x <- r[1L]; y = r[2L]
    fmatch(list(lon = c(x-0.5, x-0.5, x-0.5, x, x+0.5, x+0.5, x+0.5, x), 
                lat = c(y-0.5, y, y+0.5, y+0.5, y+0.5, y, y-0.5, y-0.5)), 
           .subset(shape, c("longitude", "latitude"))) |> na_rm()
})

table(lengths(neigh))
neigh <- data.frame(from = rep(seq_along(neigh), lengths(neigh)),
                    to   = unlist(neigh[lengths(neigh) > 0], use.names = FALSE)) 

# Remove duplicates
neigh %<>% mutate(from_to = pmin(from, to), 
                  to_from = pmax(to, from)) %>%
  unique(cols = c("from_to", "to_from")) %>% 
  subset(from != to & from > 0 & to > 0)
  
# Now we can create the adjacency matrix
if(FALSE) mapview::mapview(subset(shape, unique(c(neigh$from, neigh$to)))) +
  mapview::mapview(st_as_sf(qDF(atomic_elem(shape)), coords = c("longitude", "latitude"), crs = 4326))

if(TRUE) { # Population weighting
# Getting population centroids
WPOP <- terra::rast("/Users/sebastiankrantz/Documents/Data/WorldPop/ppp_2020_1km_Aggregated.tif")
  
# Segments
ext <- exactextractr::exact_extract(WPOP, shape, include_xy = TRUE) |> rowbind(idcol = "ID")

# Population-weighted mean of x and y within each polygon
ext %<>% 
  mutate(pop_weight = value * coverage_fraction) %>%
  group_by(ID) %>% 
  summarise(across(c(x, y), fmean, pop_weight),
            pop_weight = fsum(pop_weight)) 

if(FALSE) {  
# Centroids not in the shape
not_inside <- which(!s2::s2_within(with(ext, s2::s2_lnglat(x, y)), shape$geometry))

# Getting nearest point
nearest_points <- st_nearest_points(shape[not_inside, ],
  ext[not_inside, ] |> st_as_sf(coords = c("x", "y"), crs = 4326), pairwise = TRUE) |> 
  st_cast("POINT", group_or_split = FALSE, warn = FALSE) |> st_as_sf()

if(FALSE) ext[not_inside, ] |>
  st_as_sf(coords = c("x", "y"), crs = 4326) |>
  mapview::mapview(zcol = "pop_weight", layer.name = "Population") +
  mapview::mapview(shape[not_inside, ]) +
  mapview::mapview(nearest_points)

# Replacing 
set(ext, not_inside, c("x", "y"), mctl(st_coordinates(nearest_points)))
}

tfm(shape) <- ext |> roworder(ID) |> select(pwx = x, pwy = y, pop_weight) 
settfm(shape, pwx = pfirst(pwx, longitude), pwy = pfirst(pwy, latitude))

# Test
if(FALSE) mapview::mapview(subset(shape, unique(c(neigh$from, neigh$to)))) +
  mapview::mapview(st_as_sf(qDF(atomic_elem(shape)), coords = c("pwx", "pwy"), crs = 4326))
} # -> segments with population weighted data doesn't work too well

# Saving shape
qsave(shape, "data/grid_network/cells_shape.qs")

# Now generating segments
lon_col <- shape[["pwx"]] # "longitude"
lat_col <- shape[["pwy"]] # "latitude"

segments <- Map(
  f = function(i, j) {
    st_linestring(matrix(c(lon_col[i], lon_col[j], 
                           lat_col[i], lat_col[j]), 2))
  },
  neigh$from,
  neigh$to
)

segments <- st_sf(neigh, geometry = st_sfc(segments, crs = 4326))

if(FALSE) mapview::mapview(segments) 
# line2df(segments)
# remove <- c(8461, 8462, 8454, 8457, 8455, 8449, 8453, 8442, 8443, 
#             8435, )

# Saving segments
qsave(segments, "data/grid_network/segments.qs")


# Now scraping the routing data following Graff (2024)
routes <- data.table(from = segments$from, 
                     to = segments$to, 
                     from_lon = lon_col[segments$from],
                     from_lat = lat_col[segments$from],
                     to_lon = lon_col[segments$to],
                     to_lat = lat_col[segments$to],
                     duration = NA_real_, 
                     distance = NA_real_, 
                     geometry = list())
# Fetch Routes
for (i in seq_row(routes)) {
  cat(i, " ")
  route <- osrmRoute(ss(routes, i, c("from_lon", "from_lat")),
                     ss(routes, i, c("to_lon", "to_lat")), overview = "simplified") |>
    tryCatch(error = function(e) NULL)
  if(is.null(route)) {
    cat(sprintf("\nroute %d from %s to %s could not be calculated\n", i, routes$from_city_ascii[i], routes$to_city_ascii[i]))
    next    
  }
  if(i %% 1000 == 0) Sys.sleep(5)
  set(routes, i, 7:9, select(route, duration, distance, geometry))
}
routes <- routes |> subset(is.finite(duration)) |> st_as_sf(crs = st_crs(route))

# Adding route start and end points
add_vars(routes) <- select(line2df(routes), -L1)
# Compute distances
settfm(routes, 
  sp_distance = geodist::geodist_vec(from_lon, from_lat, to_lon, to_lat, 
       measure = "geodesic", paired = TRUE) / 1000,
  start_distance = geodist::geodist_vec(fx, fy, from_lon, from_lat, 
       measure = "geodesic", paired = TRUE) / 1000,
  end_distance = geodist::geodist_vec(tx, ty, to_lon, to_lat, 
       measure = "geodesic", paired = TRUE) / 1000
)

routes %<>% base::subset(distance > 0)

routes |> gvr("distance") |> descr()
routes %$% descr(start_distance / distance)
routes %$% descr(end_distance / distance)

settfm(routes,
  start_duration_10kmh = 60 * start_distance / 10,
  end_duration_10kmh = 60 * end_distance / 10
)

routes %$% descr(start_duration_10kmh / duration)
routes %$% descr(end_duration_10kmh / duration)

settfm(routes, 
  speed_kmh = distance / ((duration + start_duration_10kmh + end_duration_10kmh) / 60),
  route_efficiency = sp_distance / (distance + start_distance + end_distance),
  time_efficiency = sp_distance / ((duration + start_duration_10kmh + end_duration_10kmh) / 60)
)

descr(routes$route_efficiency)
descr(routes$time_efficiency)

# Saving
routes |> qsave("data/grid_network/routes_raw_pw.qs")

# Plotting
mapview::mapview(shape) + (segments |> 
  join(routes, on = c("from", "to"), how = "inner", drop = TRUE) |> 
  select(time_efficiency) |> 
  mapview::mapview())
  
  
