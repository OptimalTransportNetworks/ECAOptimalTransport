#############################################
# Adding High-Value Links to Existing Network
#############################################

# Note: This is a self-contained section due to the use of specific indices for manual adjustments. 
# It uses a previous version of the network identical to the current/generated version but where the nodes and edges are in a different (random) order. 
# The output is a spatial data frame 'add_links' with the final proposed links. These can be added to the current/generated version.

library(fastverse)
set_collapse(mask = c("manip", "helper", "special"), nthreads = 4)
fastverse_extend(qs, sf, s2, units, sfnetworks, stplanr, install = TRUE)
source("code/helpers/helpers.R")
fastverse_conflicts()

net <- qread("data/transport_network/net_discrete_final.qs")
dist_ttime_mats <- qread("data/transport_network/net_dist_ttime_mats_adj.qs")
sym_dist_mat <- (dist_ttime_mats$distances + t(dist_ttime_mats$distances)) / 2

## Plotting
plot(net)
nodes <- net |> activate("nodes") |> tidygraph::as_tibble() |> 
  mutate(.tidygraph_node_index = NULL)
edges <- net |> activate("edges") |> tidygraph::as_tibble() |> 
  mutate(.tidygraph_edge_index = NULL)

nodes_dmat <- st_distance(nodes) |> set_units("km")
diag(nodes_dmat) <- NA
nodes_dmat[upper.tri(nodes_dmat)] <- NA

# Routes to be calculated
routes_ind <- which(!is.na(nodes_dmat), arr.ind = TRUE)
nrow(routes_ind)

nodes_coord <- st_coordinates(nodes) |> qDF() |> setNames(c("lon", "lat"))

keep_routes <- !intercepted_routes(routes_ind, nodes_coord, sym_dist_mat, feasible = TRUE, alpha = 45, mr = 1/0.767, fmr = 1.5) # US: 0.843 EU: 0.767
sum(keep_routes)

add_links <- with(nodes_coord, lapply(mrtl(routes_ind[keep_routes, ]), function(r) st_linestring(cbind(lon[r], lat[r])))) |> 
              st_sfc(crs = 4326) |> st_as_sf()
add_links_df <- line2df(add_links)[-1] %*=% 1000 %>% dapply(as.integer)
edges_df <- line2df(edges)[-1] %*=% 1000 %>% dapply(as.integer)
m <- add_links_df %in% edges_df | add_links_df %in% gv(edges_df, c(3:4, 1:2))
nrow(add_links) - sum(m)
add_links <- add_links[!m, ]
rm(add_links_df, edges_df)

add_links %<>% ss(-c(265, 266, 218, 179, 171, 165, 195, 181, 117, 116, 87, 78, 
                     74, 73, 72, 59, 60, 61, 62, 258, 259, 187, 148, 147, 135, 
                     263, 268, 237, 304))

# mapview(qread("data/transport_network/segments.qs"), map.types = c(mapviewGetOption("basemaps"), "Esri.WorldStreetMap", "Esri.WorldTerrain")) + mapview(nodes) + mapview(add_links, color = "green")

add_links |> qsave("data/transport_network/add_links_network_30km_alpha45_mrEU_fmr15.qs")
