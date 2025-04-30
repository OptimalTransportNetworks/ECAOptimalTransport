library(fastverse)
set_collapse(mask = c("manip", "helper", "special"), nthreads = 4)
fastverse_extend(qs, sf, s2, units, stplanr, sfnetworks, osrm, tmap, install = TRUE)
source("code/helpers.R")
fastverse_conflicts()

####################################
# Part 1: Important Cities and Ports
####################################

countries <- wbstats::wb_cachelist$countries %$% iso3c[region_iso3c %==% "ECS"] %>% 
  c("XKO") %>% setdiff(c("GRL", "ISL", "GBR", "IRL", "NOR", "SWE", "FIN", "FRO", "IMN", "CYP"))

# From geojson.io
ECA_outline <- st_read("/Users/sebastiankrantz/Documents/ECAOptimalTransport/data/outline/layers") |> 
               st_buffer(as_units(5, "km"))

# Local GDP data: 0.5 degree, censored for density below 0.01
LGDP <- fread("/Users/sebastiankrantz/Documents/Data/LocalGDP/0_5deg/final_GDP_0_5deg_postadjust_pop_dens_0_01_adjust.csv") |> 
        subset(iso %in% countries) |>
        transform(longitude = longitude + 0.25, # Coordinate is bottom left corner
                  latitude = latitude + 0.25) 
  
LGDP %<>% subset(LGDP |> st_as_sf(coords = c("longitude", "latitude"), crs = 4326) |> 
                 st_intersects(ECA_outline, sparse = FALSE))
shape <- st_read("/Users/sebastiankrantz/Documents/Data/LocalGDP/0_5deg/shapefile") |>
         base::subset(cell_id %in% LGDP$cell_id)
shape %<>% join(LGDP, on = c("cell_id", "subcell_id"), how = "inner", drop = TRUE)
shape %<>% st_make_valid() 

mapview::mapview(shape, zcol = "pop_cell") + 
  mapview::mapview(ECA_outline, color = "red")

# Creating Adjacency grid  -----------------------------------------------

# queen = TRUE  →  edge or corner touching
neigh <- spdep::poly2nb(shape, queen = TRUE)
neigh <- data.frame(from = rep(seq_len(nrow(shape)), lengths(neigh)),
                    to   = unlist(neigh, use.names = FALSE)) 
# Remove duplicates
neigh %<>% mutate(from_to = pmin(from, to), 
                  to_from = pmax(to, from)) %>%
  unique(cols = c("from_to", "to_from")) %>% 
  roworder(from, to) %>%
  subset(from != to & from > 0 & to > 0 & rowid(from) < 9) %>% 
  roworder(to, from) %>%
  subset(rowid(to) < 9)
  
# Now we can create the adjacency matrix
mapview::mapview(subset(shape, unique(c(neigh$from, neigh$to)))) +
  mapview::mapview(st_as_sf(qDF(atomic_elem(shape)), coords = c("longitude", "latitude"), crs = 4326))

if(FALSE) {
# Getting population cetroids
WPOP <- terra::rast("/Users/sebastiankrantz/Documents/Data/WorldPop/ppp_2020_1km_Aggregated.tif")
  
# Segments
ext <- terra::extract(WPOP, shape,
                      cells   = TRUE,   # return the raster cell number
                      weights = TRUE,   # return fractional coverage
                      na.rm   = TRUE)   # drop cells where population is NA

# Append the x / y coordinates of every raster cell
xy <- terra::xyFromCell(WPOP, ext$cell)
ext$x <- xy[, 1]
ext$y <- xy[, 2]

# Population-weighted mean of x and y within each polygon
setDT(ext)                                            # data.table magic
ext[, pop_weight := ppp_2020_1km_Aggregated * weight] # weight = pop * area share

ext %<>% group_by(ID) %>% 
  summarise(across(c(x, y), fmean, pop_weight),
            pop_weight = fsum(pop_weight)) %>% 
  mutate(ID = as.integer(ID))

add_vars(shape) <- ext |> roworder(ID) |> select(pwx = x, pwy = y, pop_weight) 
settfm(shape, pwx = pfirst(pwx, longitude), pwy = pfirst(pwy, latitude))

# Test
mapview::mapview(subset(shape, unique(c(neigh$from, neigh$to)))) +
  mapview::mapview(st_as_sf(qDF(atomic_elem(shape)), coords = c("pwx", "pwy"), crs = 4326))
} # -> segments with population weighted data doesn't work too well


# Now generating segments
lon_col <- "longitude"
lat_col <- "latitude"

segments <- Map(
  f = function(i, j) {
    st_linestring(matrix(c(shape[[lon_col]][i], shape[[lon_col]][j], 
                           shape[[lat_col]][i], shape[[lat_col]][j]), 2))
  },
  neigh$from,
  neigh$to
)

segments <- st_sfc(segments, crs = 4326)
neigh_lines <- st_sf(neigh, geometry = segments)

mapview::mapview(neigh_lines) 

line2df(neigh_lines)

# Now scraping the routing data following Graff (2024)
library(osrm)
routes <- data.table(from = neigh_lines$from, 
                     to = neigh_lines$to, 
                     from_lon = shape$longitude[neigh_lines$from],
                     from_lat = shape$latitude[neigh_lines$from],
                     to_lon = shape$longitude[neigh_lines$to],
                     to_lat = shape$latitude[neigh_lines$to],
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
  route_efficiency = sp_distance / (distance + start_distance + end_distance),
  time_efficiency = sp_distance / ((duration + start_duration_10kmh + end_duration_10kmh) / 60)
)

descr(routes$route_efficiency)
descr(routes$time_efficiency)

# Saving
routes |> qsave("data/grid_network/routes_raw.qs")


