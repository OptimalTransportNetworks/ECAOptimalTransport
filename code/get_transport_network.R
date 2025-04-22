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

# Get GADM Boundaries
st_layers("/Users/sebastiankrantz/Documents/Data/GADM/gadm_410-levels.gpkg")
GADM0 <- st_read("/Users/sebastiankrantz/Documents/Data/GADM/gadm_410-levels.gpkg", layer = "ADM_0") |> 
  base::subset(GID_0 %in% countries) |> st_cast("POLYGON") |>
  base::subset(st_centroid(geom) |> st_coordinates() |> qDF() |> with(X > -12)) |>
  rmapshaper::ms_simplify(keep = 0.02) # |> st_make_valid()
setdiff(countries, GADM0$GID_0)
setdiff(GADM0$GID_0, countries)
mapview::mapview(GADM0)

# # Get cities: https://www.geonames.org/export/dump/
# cities <- fread("/Users/sebastiankrantz/Documents/Data/cities/geonames-all-cities-with-a-population-1000.csv") |>
#   janitor::clean_names() |>
#   mutate(iso3c = countrycode::countrycode(country_code, "iso2c", "iso3c"), 
#          iso3c = ifelse(country_code == "XK", "XKO", iso3c)) |>
#   transform(tstrsplit(coordinates, ",") |> set_names(c("lat", "lon")) |> lapply(as.numeric)) |>
#   subset(iso3c %in% countries & is.finite(population) & population > 2e4)

# Simplemaps is better! (more accurate populations. Some geonames very old (from 2006))
cities <- fread("/Users/sebastiankrantz/Documents/Data/cities/simplemaps_worldcities_basicv1.77/worldcities.csv") |>
          subset(iso3 %in% countries & is.finite(population)) |> rename(lng = lon)

setdiff(countries, cities$iso3)
cities <- cities |> transform(x = lon, y = lat) |> 
  st_as_sf(coords = c("x", "y"), crs = 4326) |>
  subset(fnobs(st_within(geometry, GADM0$geom)) > 0)

# mapview::mapview(cities) 
if(FALSE) { # Test: Paris Metropolitan Area
  near_paris <- which(unclass(st_distance(cities, subset(cities, city == "Paris"))) < 10000)
  sum(cities$population[near_paris])
}
# cities |> subset(population > 1e6) |> mapview::mapview() 

cities_red <- largest_within_radius(subset(cities, population >= 2e5 & lon < 83 & lat < 60), c("lon", "lat"), "population", radius_km = 30) 
cities_red <- join_within_radius(cities_red, cities, c("lon", "lat"), size = "population", radius_km = 20)
cities_red <- largest_within_radius(subset(cities_red, population >= 2.5e5), c("lon", "lat"), "population", radius_km = 50)  
# cities_red <- largest_within_radius(cities_red, c("lon", "lat"), "population", radius_km = 100)  
mapview::mapview(cities_red)

# Plotting
# cities_red |> with(plot(lon, lat, cex = sqrt(population)/1e3, pch = 16, col = rgb(0, 0, 0, alpha=0.5)))
# cities_red |> select(population, city) |> mapview::mapview()

# Remove unnavigable towns
cities_red %<>% subset(city_ascii %!in% c("Palma", "Palermo", "Catania", "Messina"))

fndistinct(atomic_elem(cities_red))

# Raw matrix of routes:
dist_ttime_mats <- split_large_dist_matrix(setRownames(select(qDF(cities_red), lon, lat), cities_red$city_ascii))

# Checks
all.equal(dist_ttime_mats$sources, dist_ttime_mats$destinations)
allv(diag(dist_ttime_mats$durations), 0)
allv(diag(dist_ttime_mats$distances), 0)
diag(cor(dist_ttime_mats$sources, select(qDF(cities_red), lon, lat)))
s2_distance(with(dist_ttime_mats$sources, s2_lnglat(lon, lat)),
            with(select(cities_red, lon, lat), s2_lnglat(lon, lat))) |> descr()
pwcor(unattrib(dist_ttime_mats$distances), unattrib(dist_ttime_mats$durations))

# Saving
cities_red |> atomic_elem() |> qDF() |> fwrite("data/transport_network/cities.csv")
dist_ttime_mats |> qsave("data/transport_network/cities_dist_ttime_mats.qs")

# Cleanup
rm(list = ls())
source("code/helpers.R")
fastverse_conflicts()


####################################
# Part 2: Transport Network
####################################

# Reading again
cities <- fread("data/transport_network/cities.csv")
dist_ttime_mats <- qread("data/transport_network/cities_dist_ttime_mats.qs")

# Distance between city centroids and route start/end points
descr(diag(st_distance(st_as_sf(cities, coords = c("lon", "lat"), crs = 4326), 
                       st_as_sf(dist_ttime_mats$sources, coords = c("lon", "lat"), crs = 4326))))
# Create route-start-point version of the dataset
cities_rsp_sf <- cities |> 
  mutate(geometry = st_as_sf(dist_ttime_mats$sources, coords = c("lon", "lat"), crs = 4326)$geometry) |>
  st_as_sf()

# Routes between all connections ------------------------------------------

cities_dmat <- dist_ttime_mats$sources |> with(s2_lnglat(lon, lat)) |> st_distance() |> set_units("km")
diag(cities_dmat) <- NA
cities_dmat[upper.tri(cities_dmat)] <- NA
cities_dmat[cities_dmat > as_units(2000, "km")] <- NA

if(FALSE) {
# Create a leaflet map
library(leaflet)
leaflet(data = cities_rsp_sf) %>%
  addTiles() %>%  # Add base map tiles
  addCircleMarkers(
    ~lon, ~lat,  # Coordinates from df
    popup = ~paste("Population:", population),  # Popup with population info
    radius = 8,  # Marker size
    color = "blue",
    fillOpacity = 0.6
  ) %>%
  addMeasure(
    position = "topright",  # Position of the measure control
    primaryLengthUnit = "meters",  # Measure distances in meters
    primaryAreaUnit = "sqmeters",  # Measure areas in square meters
    activeColor = "#FF0000",  # Color when measurement is active
    completedColor = "#00FF00"  # Color when measurement is completed
  )
}

# Routes to be calculated
routes_ind <- which(!is.na(cities_dmat), arr.ind = TRUE)
nrow(routes_ind)

# Determining Ideal Hypothetical (US-Grade) Network 
# See: https://christopherwolfram.com/projects/efficiency-of-road-networks/

# EU Route Efficiency  
keep_routes <- !intercepted_routes(routes_ind, dist_ttime_mats$sources, NULL, alpha = 45, mr = 1/0.767) 
sum(keep_routes)

# Plot ideal Network
with(cities, {
  oldpar <- par(mar = c(0,0,0,0))
  plot(lon, lat, cex = sqrt(population)/1e3, pch = 16, col = rgb(0, 0, 0, alpha=0.5), 
       axes = FALSE, xlab = NA, ylab = NA, asp=1)
  for (r in mrtl(routes_ind[keep_routes, ])) { # Comment loop to get LHS of Figure 14
    lines(lon[r], lat[r])
  }
  par(oldpar)
})


dev.copy(pdf, "figures/trans_ECA_network_EU_45deg.pdf", width = 10, height = 3.5)
dev.off()


# Fetching all Routes and Simplifying -------------------------------------------

# -> Better skip and load finished segments below

routes <- data.table(from_city_ascii = cities$city_ascii[routes_ind[, 1]], 
                     to_city_ascii = cities$city_ascii[routes_ind[, 2]], 
                     duration = NA_real_, 
                     distance = NA_real_, 
                     geometry = list())
# Fetch Routes
for (i in seq_row(routes_ind)) {
  cat(i, " ")
  route <- osrmRoute(ss(cities, routes_ind[i, 1L], c("lon", "lat")),
                     ss(cities, routes_ind[i, 2L], c("lon", "lat")), overview = "full") |>
        tryCatch(error = function(e) NULL)
  if(is.null(route)) {
    cat(sprintf("\nroute %d from %s to %s could not be calculated\n", i, routes$from_city_ascii[i], routes$to_city_ascii[i]))
    next    
  }
  if(i %% 1000 == 0) Sys.sleep(10)
  set(routes, i, 3:5, select(route, duration, distance, geometry))
}
routes <- routes |> subset(is.finite(duration)) |> st_as_sf(crs = st_crs(route))

# Saving
routes |> qsave("data/transport_network/routes_raw.qs")

# Adding Gravity to Routes (https://en.wikipedia.org/wiki/Newton%27s_law_of_universal_gravitation)
# routes <- qread("data/transport_network/routes_raw.qs")
dmat <- st_distance(cities_rsp_sf)
diag(dmat) <- NA
frange(dmat) 

colnames(dmat) <- rownames(dmat) <- cities_rsp_sf$city_ascii
ddf <- pivot(qDF(dmat, "from"), "from", names = list("to", "sp_dist"), na.rm = TRUE) |> 
  join(select(cities, from = city_ascii, from_pop = population)) |> 
  join(select(cities, to = city_ascii, to_pop = population)) |> 
  mutate(gravity = as.double(from_pop) * to_pop / sp_dist / 1e6) # Figure in TEU for port cities??

routes <- routes |> 
  join(ddf, on = c(from_city_ascii = "from", to_city_ascii = "to"), drop = "x") |>
  mutate(gravity_rd = as.double(from_pop) * to_pop / distance / 1e6, 
         gravity_dur = as.double(from_pop) * to_pop / duration / 1e6)

# Intersecting Routes
segments <- lapply(gsplit(g = sample.int(10, nrow(routes), TRUE)), function(ind) {
     overline2(ss(routes, ind) |> mutate(passes = 1L), 
        attrib = c("passes", "gravity", "gravity_rd", "gravity_dur"))
})
rm(routes); gc()
segments %<>% rowbind() %>% ss(!is_linepoint(.)) %>% st_make_valid()
segments %<>% overline2(attrib = c("passes", "gravity", "gravity_rd", "gravity_dur"))
gc()

# Saving
segments |> qsave("data/transport_network/segments.qs")


# Loading Segments --------------------------------------------------------------------

segments <- qread("data/transport_network/segments.qs") 

# First Round of subdivision
segments <- rmapshaper::ms_simplify(segments, keep = 0.2, snap_interval = deg_m(500)) |> 
  subset(vlengths(geometry) >= 4) # |> st_cast("LINESTRING") |> st_make_valid() |> st_cast("LINESTRING")
segments <- overline2(segments, attrib = c("passes", "gravity", "gravity_rd", "gravity_dur"))
if(any(is_linepoint(segments))) segments %<>% ss(!is_linepoint(.))
# plot(segments)
# mapview::mapview(segments) + mapview::mapview(rnet_get_nodes(segments))

# -> Feel free to skip and load final network below

# Creating Network
net <- as_sfnetwork(segments, directed = FALSE)
# plot(net)
summ_fun <- function(fun) list(passes = fun, gravity = fun, gravity_rd = fun, gravity_dur = fun, "ignore")
filter_smooth <- function(net) {
  net |> 
    tidygraph::activate("edges") |> 
    dplyr::filter(!tidygraph::edge_is_multiple()) |> 
    dplyr::filter(!tidygraph::edge_is_loop()) |> 
    tidygraph::convert(to_spatial_smooth, 
                       protect = cities_rsp_sf$geometry,
                       summarise_attributes = summ_fun("mean"))
}


## Filtering, smoothing and spatial subdivision 
net <- filter_smooth(net)
net <- filter_smooth(net)
# plot(net)
# mapview(st_geometry(net, "edges")) + mapview(st_geometry(net, "nodes"))

# Saving Smoothed Version (pre contraction)
net |> qsave("data/transport_network/net_smoothed.qs")
net <- qread("data/transport_network/net_smoothed.qs")

## Contracting network: Manually 
segments <- net |> activate("edges") |> tidygraph::as_tibble() |> mutate(.tidygraph_edge_index = NULL)
nodes <- nodes_max_passes(segments, attrib = c("passes", "gravity", "gravity_rd", "gravity_dur"))
nodes_clustered <- cluster_nodes_by_cities(nodes, select(cities_rsp_sf, -lon, -lat), 
                                           city_radius_km = 33, 
                                           cluster_radius_km = 28, 
                                           algo.dbscan = FALSE,
                                           weight = "gravity_rd")
segments_contracted <- contract_segments(segments, nodes_clustered, attrib = c("passes", "gravity", "gravity_rd", "gravity_dur"))
# plot(segments_contracted)
# mapview(segments_contracted) + mapview(rnet_get_nodes(segments_contracted))

net <- as_sfnetwork(segments_contracted, directed = FALSE)

# Filtering, smoothing and spatial subdivision 
net <- filter_smooth(net)
net <- filter_smooth(net)

# plot(net)
# mapview(net |> activate("edges") |> tidygraph::as_tibble() |> mutate(.tidygraph_edge_index = NULL)) + 
#   mapview(st_as_sf(st_geometry(net, "nodes")))

## Saving 
net |> qsave("data/transport_network/net_discrete_final.qs")

## Loading Final Network -----------------------------------------

net <- qread("data/transport_network/net_discrete_final.qs")

## Plotting
plot(net)
nodes <- net |> activate("nodes") |> tidygraph::as_tibble() |> 
  mutate(.tidygraph_node_index = NULL)
edges <- net |> activate("edges") |> tidygraph::as_tibble() |> 
  mutate(.tidygraph_edge_index = NULL, log10_gravity = log10(gravity))

nodes$key_city <- rowSums(st_distance(nodes, cities_rsp_sf) < as_units(20, "km")) > 0
sum(nodes$key_city)
descr(edges$gravity_rd)

# Needed throughout 
nodes_coord <- st_coordinates(nodes) |> qDF() |> setNames(c("lon", "lat"))

tmap_mode("plot")
tmap_options(raster.max_cells = 1e7)


# First a plot of just the routes
pdf("figures/trans_ECA_network_actual.pdf", width = 10, height = 4.2)
tm_basemap("Esri.WorldGrayCanvas", zoom = 5) +
  tm_shape(segments) + tm_lines(col = "black") +
  tm_shape(subset(nodes, key_city)) + tm_dots(size = 0.2, fill = "orange2") +
  tm_layout(frame = FALSE)
dev.off()

# Now the plot with the discretized representation
# <Figure 12: RHS>
pdf("figures/trans_ECA_network_actual_discretized_gravity_plus_orig.pdf", width = 10, height = 4.2)
tm_basemap("Esri.WorldGrayCanvas", zoom = 5) +
  tm_shape(segments) + tm_lines(col = "black") +
  tm_shape(mutate(edges, gravity_rd = gravity_rd / 1e4)) +
  tm_lines(col = "gravity_rd",
           col.scale = tm_scale_intervals(values = "inferno", breaks = c(0, 1, 5, 25, 100, Inf)*2),
           col.legend = tm_legend("Sum of Gravity", position = c("left", "top"), frame = FALSE, 
                                  text.size = 1, title.size = 1.2, 
                                  title.padding = c(-0.3, 0, -0.5, 0), 
                                  item.space = 0), lwd = 1.5) +
  tm_shape(subset(nodes, key_city)) + tm_dots(size = 0.2) +
  tm_shape(subset(nodes, !key_city)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE)
dev.off()


# Recalculating Routes Matrix for Network ----------------------------------------------------------

# -> Feel free to skip and load matrix below

nrow(nodes_coord)
all.equal(unattrib(nodes_coord), mctl(st_coordinates(nodes)))
# This can take a few minutes: now generating distance matrix of all nodes
dist_ttime_mats <- split_large_dist_matrix(nodes_coord, verbose = TRUE)
# Checks
all.equal(dist_ttime_mats$sources, dist_ttime_mats$destinations)
allv(diag(dist_ttime_mats$durations), 0)
allv(diag(dist_ttime_mats$distances), 0)
diag(cor(dist_ttime_mats$sources, nodes_coord))
s2_distance(with(dist_ttime_mats$sources, s2_lnglat(lon, lat)),
            with(nodes_coord, s2_lnglat(lon, lat))) |> descr()
pwcor(unattrib(dist_ttime_mats$distances), unattrib(dist_ttime_mats$durations))
# Now finding places that are on islands (e.g. Zanzibar City): should not exist here
if(any(which(fnobs(dist_ttime_mats$durations) < 20))) stop("Found Islands")

dist_ttime_mats |> qsave("data/transport_network/net_dist_ttime_mats_adj.qs")

# Loading again -----------------------------------

dist_ttime_mats <- qread("data/transport_network/net_dist_ttime_mats_adj.qs")

# Check
all.equal(st_as_sf(net, "nodes")$geometry, nodes$geometry)

# Making symmetric
sym_dist_mat <- (dist_ttime_mats$distances + t(dist_ttime_mats$distances)) / 2
sym_time_mat <- (dist_ttime_mats$durations + t(dist_ttime_mats$durations)) / 2

# Add average distance and travel time to edges
edges_ind <- edges |> qDF() |> select(from, to) |> qM()
edges$sp_distance <- st_length(edges)
edges$distance <- sym_dist_mat[edges_ind]
edges$duration <- sym_time_mat[edges_ind]

descr(with(edges, (distance / 1000) / (duration / 60)))

# Loading, Integrating and Plotting Additional Connections ---------------------------------------------------

# -> To generate the 'add_links' file, execute '7.1_add_links.R' in a clean R session

add_links <- qread("data/transport_network/add_links_network_30km_alpha45_mrEU_fmr15.qs")
add_links_df <- line2points(add_links)
dmat <- st_distance(nodes$geometry, add_links_df$geometry)
add_links_df$node <- dapply(dmat, which.min)
add_links_df$geometry <- nodes$geometry[add_links_df$node]
add_links <- add_links_df |> group_by(id) |> 
  summarise(from = ffirst(node),
            to = flast(node),
            geometry = st_combine(geometry)) |> st_cast("LINESTRING")
# Checks
all(line2df(add_links) %>% select(fx, fy) %in% qDF(st_coordinates(nodes)))
all(line2df(add_links) %>% select(tx, ty) %in% qDF(st_coordinates(nodes)))
all(select(line2df(add_links), fx:ty) %!in% select(line2df(edges), fx:ty))
rm(dmat, add_links_df)

tmap_mode("plot")

# Same as Figure 15 but with discrete edges (Figure 15 with real roads is obtained below)
pdf("figures/trans_ECA_network_actual_discretized_gravity_new_roads_Esri.WorldStreetMap.pdf", width = 10, height = 4.2)
tm_basemap("Esri.WorldTopoMap", zoom = 5) +
  tm_shape(segments) +
  tm_lines(col = "gravity_rd",
           col.scale = tm_scale_intervals(values = "inferno", breaks = c(0, 1, 5, 25, 100, Inf)),
           col.legend = tm_legend("Sum of Gravity", position = c("right", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) +
  tm_shape(add_links) + tm_lines(col = "green4", lwd = 1) + 
  tm_shape(subset(nodes, key_city)) + tm_dots(size = 0.2) +
  tm_shape(subset(nodes, !key_city)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE)
dev.off()



# Now Adding Populations ----------------------------------------------------------

# Distance Matrix
dmat <- nodes |> st_distance(st_as_sf(cities, coords = c("lon", "lat"), crs = 4326)) 
# Remove cities part of port cities
used <- dapply(dmat[nodes$key_city, ] < as_units(30, "km"), any)
table(used)
dmat <- dmat[!nodes$key_city, !used]
dim(dmat)
dmat[dmat >= as_units(30, "km")] <- NA
col <- numeric(nrow(dmat))
m <- dapply(dmat, function(x) if (allNA(x)) col else replace(col, which.min(x), 1)) # Assigning each city to closest node
nodes$population <- NA
nodes$population[!nodes$key_city] <- m %*% cities$population[!used] # Summing populations of closest cities
nodes$city_ascii <- NA_character_
nodes$city_ascii[!nodes$key_city] <- qDT(cities) |> 
  extract(!used, paste(city_ascii, iso3, sep = " - ")) |> 
  extract(dapply(m %r*% cities$population[!used], 
                 function(x) if(any(x > 0.1)) which.max(x) else NA_integer_, MARGIN = 1))
# Now adding port cities (closest matches)
ind <- dapply(st_distance(cities_rsp_sf, nodes), function(x) if (any(x < 20e3)) which.min(x) else NA_integer_)
nodes$population[!is.na(ind)] <- cities_rsp_sf$population[na_rm(ind)]
nodes$city_ascii[!is.na(ind)] <- cities_rsp_sf$city_ascii[na_rm(ind)]
# Cleanup
rm(col, m, ind, used, dmat)

# Ratios
sum(nodes$population) / sum(cities_rsp_sf$population)
sum(nodes$population) / sum(cities$population)
(sum(nodes$population > 0) - nrow(cities_rsp_sf)) / (nrow(nodes) - nrow(cities_rsp_sf))

# Saving all necessary objects in an RData file ---------------------------------------------------

# First adding info to edges
tmp <- net |> st_as_sf("edges") |> atomic_elem() |> qDT() |> 
  join(atomic_elem(edges), on = c("from", "to", "passes"), overid = 2L, drop.dup.cols = "x") |> 
  select(-from, -to, -.tidygraph_edge_index)

net %<>% activate("edges") %>% dplyr::mutate(tmp)
rm(tmp)

save(nodes, edges, edges_ind, nodes_coord, net, add_links,  
     cities, cities_rsp_sf, 
     dist_ttime_mats, sym_dist_mat, sym_time_mat, 
     file = "data/transport_network/trans_ECA_network.RData")


##########################################################
# Fetching Simplified Routes (Edges) for Visual Inspection
##########################################################

load("data/transport_network/trans_ECA_network.RData")

# This is the previous Plot
tm_basemap("Esri.WorldGrayCanvas", zoom = 5) +
  tm_shape(mutate(edges, gravity_rd = gravity_rd / 1e4)) +
  tm_lines(col = "gravity_rd",
           col.scale = tm_scale_intervals(values = "inferno", breaks = c(0, 1, 5, 25, 100, Inf)*2),
           col.legend = tm_legend("Sum of Gravity", position = c("left", "top"), frame = FALSE, 
                                  text.size = 1, title.size = 1.2, 
                                  title.padding = c(-0.3, 0, -0.5, 0), 
                                  item.space = 0), lwd = 1.5) +
  tm_shape(add_links) + tm_lines(col = "green4", lwd = 1) + 
  tm_shape(subset(nodes, key_city)) + tm_dots(size = 0.2) +
  tm_shape(subset(nodes, !key_city)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE) 

# dev.copy(pdf, "figures/trans_ECA_network_actual_discretized_gravity_new_roads.pdf", 
#          width = 10, height = 10)
# dev.off()

# Fetching Simplified Routes

# -> Feel free to skip and load result below

edges_real <- edges |> qDT() |> 
  transform(geometry = list(NULL), distance = NA_real_, duration = NA_real_)
nodes_coord <- st_coordinates(nodes) |> qDF() |> setNames(c("lon", "lat"))

# Fetch Routes
for (i in seq_row(edges_ind)) {
  cat(i, " ")
  # Sys.sleep(0.1)
  route <- osrmRoute(ss(nodes_coord, edges_ind[i, 1]),
                     ss(nodes_coord, edges_ind[i, 2]), overview = "simplified") |>
    tryCatch(error = function(e) NULL)
  if(is.null(route)) {
    cat(sprintf("\nroute %d from %s to %s could not be calculated", i, nodes$city_ascii[edges$from[i]], nodes$city_ascii[edges$to[i]]))
    next    
  }
  set(edges_real, i, c("duration", "distance", "geometry"), 
      select(route, duration, distance, geometry))
}
edges_real <- edges_real |> st_as_sf(crs = st_crs(route)) |> st_make_valid()
edges_real |> qsave("data/transport_network/edges_real.qs")
rm(route, i)

edges_real %<>% rmapshaper::ms_simplify(keep = 0.1)
edges_real %<>% st_make_valid()
edges_real |> qsave("data/transport_network/edges_real_simplified.qs")

# Update info
edges_real <- qread("data/transport_network/edges_real.qs")
edges_real %<>% join(atomic_elem(edges), on = c("from", "to"), drop = "x", overid = 2)

# Draw the updated plot

tmap_options(raster.max_cells = 1e8)
pdf("figures/trans_ECA_network_actual_discretized_gravity_new_roads_real_edges_OTM.pdf", width = 10, height = 4.2)
tm_basemap("OpenTopoMap", zoom = 5) + # CartoDB.Positron
  tm_shape(mutate(edges_real, gravity_rd = gravity_rd / 1e4)) +
  tm_lines(col = "gravity_rd",
           col.scale = tm_scale_intervals(values = "inferno", breaks = c(0, 1, 5, 25, 100, Inf)*2),
           col.legend = tm_legend("Sum of Gravity", position = c("left", "top"), frame = FALSE, 
                                  text.size = 1, title.size = 1.2, 
                                  title.padding = c(-0.3, 0, -0.5, 0), 
                                  item.space = 0), lwd = 0.85) +
  tm_shape(add_links) + tm_lines(col = "green4", lwd = 0.7) + # limegreen # , lty = "twodash"
  tm_shape(subset(nodes, key_city)) + tm_dots(size = 0.12, fill = "dodgerblue4") +
  tm_shape(subset(nodes, !key_city & population > 0)) + tm_dots(size = 0.08, fill = "grey20") +
  tm_shape(subset(nodes, !key_city & population <= 0)) + tm_dots(size = 0.08, fill = "grey70") +
  tm_layout(frame = FALSE) #, inner.margins = c(0.1, 0.1, 0.1, 0.1))
dev.off()


