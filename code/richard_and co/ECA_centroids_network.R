library(fastverse)
set_collapse(mask = c("manip", "helper", "special"), nthreads = 4)
fastverse_extend(qs, sf, s2, units, stplanr, sfnetworks, osrm, tmap, mapview, install = TRUE)
source("code/fs_build/helpers.R")
fastverse_conflicts()

ECA_centroids <- fread("data/ECA_centroids.csv") |> 
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

dmat <- st_distance(ECA_centroids)
fquantile(dmat)/1000
diag(dmat) <- NA
dmat[upper.tri(dmat)] <- NA
dmat[dmat > as_units(1200, "km")] <- NA

routes_ind <- which(is.finite(dmat), arr.ind = TRUE)
nrow(routes_ind)

routes <- data.table(from_shapeISO = ECA_centroids$shapeISO[routes_ind[, 1L]], 
                     to_shapeISO = ECA_centroids$shapeISO[routes_ind[, 2L]], 
                     duration = NA_real_, 
                     distance = NA_real_, 
                     geometry = list())

coord <- st_coordinates(ECA_centroids)

# Fetch Routes
for (i in whichNA(routes$duration)) {
  cat(i, " ")
  route <- osrmRoute(coord[routes_ind[i, 1L], ],
                     coord[routes_ind[i, 2L], ], overview = "full") |>
    tryCatch(error = function(e) NULL)
  if(is.null(route)) {
    cat(sprintf("\nroute %d from %s to %s could not be calculated\n", i, routes$from_shapeISO[i], routes$to_shapeISO[i]))
    next    
  }
  if(i %% 5000 == 0) Sys.sleep(60)
  set(routes, i, 3:5, fselect(route, duration, distance, geometry))
}

routes <- routes |> fsubset(is.finite(duration)) |> st_as_sf(crs = st_crs(route))

# Saving
routes |> qsave("data/ECA_centroids_network/routes_raw.qs")

# Adding Gravity to Routes (https://en.wikipedia.org/wiki/Newton%27s_law_of_universal_gravitation)
routes <- qread("data/ECA_centroids_network/routes_raw.qs")

# Intersecting Routes
segments <- lapply(gsplit(g = sample.int(10, nrow(routes), TRUE)), function(ind) {
  overline2(ss(routes, ind) |> mutate(passes = 1L), attrib = "passes")
})
rm(routes); gc()
segments %<>% rowbind() %>% ss(!is_linepoint(.)) %>% st_make_valid()
segments %<>% overline2(attrib = "passes")
gc()


# Saving
segments |> qsave("data/ECA_centroids_network/segments.qs")


# Loading Segments --------------------------------------------------------------------

segments <- qread("data/ECA_centroids_network/segments.qs") 

# First Round of subdivision
segments <- overline2(segments, attrib = "passes")
if(any(is_linepoint(segments))) segments %<>% ss(!is_linepoint(.))
segments %<>% st_make_valid()
# plot(segments)
# mapview::mapview(segments) + mapview::mapview(rnet_get_nodes(segments))

# -> Feel free to skip and load final network below

# Creating Network
net <- as_sfnetwork(segments, directed = FALSE)
# plot(net)
summ_fun <- function(fun) list(passes = fun, "ignore")
filter_smooth <- function(net) {
  net |> 
    tidygraph::activate("edges") |> 
    dplyr::filter(!tidygraph::edge_is_multiple()) |> 
    dplyr::filter(!tidygraph::edge_is_loop()) |> 
    tidygraph::convert(to_spatial_smooth, 
                       protect = ECA_centroids$geometry,
                       summarise_attributes = summ_fun("mean"))
}


## Filtering, smoothing and spatial subdivision 
net <- filter_smooth(net)
net <- tidygraph::convert(net, to_spatial_subdivision)
net <- filter_smooth(net)
# plot(net)
# mapview(st_geometry(net, "edges")) + mapview(st_geometry(net, "nodes"))

## Contracting network: Manually 
segments <- net |> activate("edges") |> tidygraph::as_tibble() |> mutate(.tidygraph_edge_index = NULL)
nodes <- nodes_max_passes(segments, attrib = "passes")
nodes_clustered <- cluster_nodes_by_cities(nodes, ECA_centroids, 
                                           city_radius_km = 33, 
                                           cluster_radius_km = 28, 
                                           algo.dbscan = FALSE,
                                           weight = "passes")
segments_contracted <- contract_segments(segments, nodes_clustered, attrib = "passes")
# plot(segments_contracted)
# mapview(segments_contracted) + mapview(rnet_get_nodes(segments_contracted))

net <- as_sfnetwork(segments_contracted, directed = FALSE)

# Filtering, smoothing and spatial subdivision 
net <- filter_smooth(net)
net <- tidygraph::convert(net, to_spatial_subdivision)
net <- filter_smooth(net)

# plot(net)
mapview(net |> activate("edges") |> tidygraph::as_tibble() |> mutate(.tidygraph_edge_index = NULL)) +
  mapview(st_as_sf(st_geometry(net, "nodes"))) + mapview(ECA_centroids)

## Saving 
net |> qsave("data/ECA_centroids_network/net_discrete_final.qs")

## Loading Final Network -----------------------------------------

net <- qread("data/ECA_centroids_network/net_discrete_final.qs")

## Plotting
plot(net)
nodes <- net |> activate("nodes") |> tidygraph::as_tibble() |> 
  mutate(.tidygraph_node_index = NULL)
edges <- net |> activate("edges") |> tidygraph::as_tibble() |> 
  mutate(.tidygraph_edge_index = NULL)

nodes$key_city <- rowSums(st_distance(nodes, ECA_centroids) < as_units(20, "km")) > 0
sum(nodes$key_city)

# Needed throughout 
nodes_coord <- st_coordinates(nodes) |> qDF() |> setNames(c("lon", "lat"))

tmap_mode("plot")
tmap_options(raster.max_cells = 1e6)


# First a plot of just the routes
pdf("figures/ECA_centroids_network_actual.pdf", width = 10, height = 4.2)
tm_basemap("Esri.WorldGrayCanvas", zoom = 5) +
  tm_shape(segments) + tm_lines(col = "black") +
  tm_shape(subset(nodes, key_city)) + tm_dots(size = 0.2, fill = "orange2") +
  tm_layout(frame = FALSE)
dev.off()

# Now the plot with the discretized representation
# <Figure 12: RHS>
pdf("figures/ECA_centroids_network_actual_discretized_plus_orig.pdf", width = 10, height = 4.2)
tm_basemap("Esri.WorldGrayCanvas", zoom = 5) +
  tm_shape(segments) + tm_lines(col = "black") +
  tm_shape(edges) +
  tm_lines(col = "passes",
           col.scale = tm_scale_intervals(values = "inferno", breaks = c(0, 1, 5, 25, 100, Inf)*2),
           col.legend = tm_legend("Passes", position = c("left", "top"), frame = FALSE, 
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

dist_ttime_mats |> qsave("data/ECA_centroids_network/net_dist_ttime_mats_adj.qs")

# Loading again -----------------------------------

dist_ttime_mats <- qread("data/ECA_centroids_network/net_dist_ttime_mats_adj.qs")

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

add_links <- qread("data/ECA_centroids_network/add_links_network_alpha45_mrEU_fmr15.qs")
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
# tmap_options(raster.max_cells = 1e8)

# # Same as Figure 15 but with discrete edges (Figure 15 with real roads is obtained below)
# pdf("figures/ECA_centroids_network_actual_discretized_new_roads_EWTM.pdf", width = 17, height = 7)
# tm_basemap("Esri.WorldTopoMap", zoom = 6) +
#   tm_shape(segments) +
#   tm_lines(col = "passes",
#            col.scale = tm_scale_intervals(values = "inferno", style = "quantile"),
#            col.legend = tm_legend_hide()) +
#   tm_shape(add_links) + tm_lines(col = "green4", lwd = 1) + 
#   tm_shape(subset(nodes, key_city)) + tm_dots(size = 0.2) +
#   tm_shape(subset(nodes, !key_city)) + tm_dots(size = 0.1, fill = "grey70") +
#   tm_layout(frame = FALSE)
# dev.off()



# Saving all necessary objects in an RData file ---------------------------------------------------

# First adding info to edges
tmp <- net |> st_as_sf("edges") |> atomic_elem() |> qDT() |> 
  join(atomic_elem(edges), on = c("from", "to", "passes"), overid = 2L, drop.dup.cols = "x") |> 
  select(-from, -to, -.tidygraph_edge_index)

net %<>% activate("edges") %>% dplyr::mutate(tmp)
rm(tmp)

save(nodes, edges, edges_ind, nodes_coord, net, add_links,  
     ECA_centroids, dist_ttime_mats, sym_dist_mat, sym_time_mat, 
     file = "data/ECA_centroids_network/ECA_centroids_network.RData")


##########################################################
# Fetching Simplified Routes (Edges) for Visual Inspection
##########################################################

load("data/ECA_centroids_network/ECA_centroids_network.RData")

# This is the previous Plot
tm_basemap("Esri.WorldGrayCanvas", zoom = 5) +
  tm_shape(edges) +
  tm_lines(col = "passes",
           col.scale = tm_scale_intervals(values = "inferno", style = "quantile"),
           col.legend = tm_legend_hide()) +
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
    cat(sprintf("\nroute %d from %s to %s could not be calculated", i, nodes$shapeISO[edges$from[i]], nodes$shapeISO[edges$to[i]]))
    next    
  }
  set(edges_real, i, c("duration", "distance", "geometry"), 
      select(route, duration, distance, geometry))
}

edges_real <- edges_real |> st_as_sf(crs = st_crs(route)) |> st_make_valid()
edges_real |> qsave("data/ECA_centroids_network/edges_real.qs")
rm(route, i)

edges_real %<>% rmapshaper::ms_simplify(keep = 0.1)
edges_real %<>% st_make_valid()
edges_real |> qsave("data/ECA_centroids_network/edges_real_simplified.qs")

# Update info
edges_real <- qread("data/ECA_centroids_network/edges_real.qs")
edges_real %<>% join(atomic_elem(edges), on = c("from", "to"), drop = "x", overid = 2)

# Draw the updated plot

tmap_options(raster.max_cells = 1e8)
pdf("figures/ECA_centroids_network_actual_discretized_new_roads_EWTM.pdf", width = 17, height = 7)
tm_basemap("Esri.WorldTopoMap", zoom = 6) +
  tm_shape(edges_real) +
  tm_lines(col = "passes",
           col.scale = tm_scale_intervals(values = "inferno", style = "quantile"),
           col.legend = tm_legend_hide()) +
  tm_shape(add_links) + tm_lines(col = "green4", lwd = 0.7) + # limegreen # , lty = "twodash"
  tm_shape(subset(nodes, key_city)) + tm_dots(size = 0.2) +
  tm_shape(subset(nodes, !key_city)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE) #, inner.margins = c(0.1, 0.1, 0.1, 0.1))
dev.off()




