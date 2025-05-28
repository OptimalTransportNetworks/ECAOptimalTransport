###############################################################################
# Classify Nodes for GE Network Simulation
###############################################################################

library(fastverse)
set_collapse(mask = c("manip", "helper", "special"))
fastverse_extend(qs, sf, sfnetworks, tmap, install = TRUE)
fastverse_conflicts()

load("data/transport_network/trans_ECA_network.RData")
edges_real <- qread("data/transport_network/edges_real_simplified.qs")


# Plot high gravity roads
tm_basemap("CartoDB.Positron", zoom = 5) +
  tm_shape(subset(edges_real, gravity_rd >= 50e3)) +
  tm_lines(col = "gravity_rd",
           col.scale = tm_scale_intervals(values = "brewer.yl_or_rd", breaks = c(0, 1, 5, 25, 100, Inf)*1e4),
           col.legend = tm_legend("Sum of Gravity", position = c("right", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) +
  tm_shape(subset(edges_real, gravity_rd < 50e3)) +
  tm_lines(col = "grey50", lwd = 2) +
  tm_shape(subset(nodes, key_city)) + tm_dots(size = 0.15) +
  tm_shape(subset(nodes, !key_city & population > 0)) + tm_dots(size = 0.1, fill = "grey20") +
  tm_shape(subset(nodes, !key_city & population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE) 

# 25 largest cities
large_cities <- nodes %$% which((population > 2.5e6 & !city_ascii %ilike% "bursa|konya|adana|kayseri|diyarbakir") | city_ascii %ilike% "copenhagen|lisbon|amsterdam|warsaw")
length(large_cities)
(largest <- nodes$city_ascii[large_cities])

# Fastest Routes between them
# igraph::all_shortest_paths(net, large_cities[1], large_cities)
large_city_paths <- lapply(large_cities, function(i) 
  st_network_paths(net, from = i, to = large_cities, weights = "duration", mode = "all")) |>
  rowbind()

large_city_paths <- list(nodes = unique(unlist(large_city_paths$node_paths)), 
                         edges = unique(unlist(large_city_paths$edge_paths)))

tm_basemap("CartoDB.Positron", zoom = 5) +
  tm_shape(subset(edges_real, large_city_paths$edges)) +
  tm_lines(col = "gravity_rd",
           col.scale = tm_scale_intervals(values = "inferno", breaks = c(0, 1, 5, 25, 100, Inf)*1e4),
           col.legend = tm_legend("Sum of Gravity", position = c("right", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) +
  tm_shape(subset(edges_real, -large_city_paths$edges)) +
  tm_lines(col = "grey70", lwd = 2) +
  tm_shape(subset(nodes, large_city_paths$nodes)) + tm_dots(size = 0.2) +
  tm_shape(subset(nodes, -large_city_paths$nodes)) + tm_dots(size = 0.1, fill = "grey50") +
  tm_shape(subset(nodes, large_cities)) + tm_dots(size = 0.2, fill = "red") +
  tm_layout(frame = FALSE) 


# Plotting -----------------------------------------------------------------------

# Classification
settfm(edges, speed_kmh = (distance / 1000) / (duration / 60))
descr(edges$speed_kmh)
tfm(edges_real) <- atomic_elem(edges) 

settfm(nodes, major_key_city = replace_na((population > 2.5e6 & !city_ascii %ilike% "bursa|konya|adana|kayseri|diyarbakir") | city_ascii %ilike% "copenhagen|lisbon|amsterdam|warsaw"))
sum(nodes$major_key_city)
descr(nodes$population[-large_cities])

settfm(nodes, product = nif(major_key_city, NA_integer_, # Heterogeneous products
                            population > 1e6, 4L,   # Metropole
                            population > 2.5e5, 3L, # Large City
                            population > 1e5, 2L,   # Medium-Sized City
                            default = 1L))          # Town/Node
table(nodes$product, na.exclude = FALSE)
setv(nodes$product, whichNA(nodes$product), seq_along(largest) + 4L)

# Need to write this to ensure product classification is available for GE simulation !!!!
nodes_param <- fread("data/transport_network/graph_nodes.csv")
nodes_param |> select(lon, lat) |> qM() |> subtract(st_coordinates(nodes)) |> descr() |> print(digits = 7)
nodes_param$product <- nodes$product
nodes_param$iso3c <- nodes$iso3c
# nodes_param |> atomic_elem() |> qDT() |> fwrite("data/transport_network/graph_nodes.csv")
rm(nodes_param)
attr(nodes$product, "levels") <- c("Small Town", "City > 100K", "City > 250K", "City > 1M", paste("Major City", seq_along(largest)))
class(nodes$product) <- "factor"

# Plotting
# <Figure 41: LHS> (Use nname <- "all_routes" above to generate the RHS)
pdf("figures/GE/trans_ECA_network_29_products.pdf", width = 10, height = 4.2)
tm_basemap("CartoDB.Positron", zoom = 5) +
  tm_shape(subset(edges_real, speed_kmh > 10)) + # Use edges_real to generate <Figure A21>
  tm_lines(col = "speed_kmh", 
           col.legend = tm_legend("Speed", position = c("right", "bottom"), stack = "h", 
                                  frame = FALSE, text.size = 1.3, title.size = 1.6),
           col.scale = tm_scale_continuous(6, values = "turbo"), # 7, 
           lwd = 2) + 
  tm_shape(droplevels(mutate(nodes, product = fifelse(unclass(product) > 4L, NA, product)))) + 
  tm_dots(fill = "product", size = 0.25, 
          fill.scale = tm_scale_categorical(values = "turbo", value.na = "purple3", label.na = "Own Product"),
          fill.legend = tm_legend("Product", position = c("right", "bottom"), # stack = "h", 
                                  frame = FALSE, text.size = 1.3, title.size = 1.6)) +
  tm_shape(subset(nodes, unclass(product) > 5L)) + tm_dots(size = 0.5, fill = "purple3") +
  tm_layout(frame = FALSE)
dev.off()         


# Projects -----------------------------------------------------------------------------------------------------

if(!all.equal(edges$gravity_rd, edges_real$gravity_rd)) stop("Mismatch!")
tfm(edges_real) <- atomic_elem(edges)
edges_real$id <- seq_row(edges_real)

mapview(edges_real, map.types = c("Esri.WorldStreetMap", mapviewGetOption("basemaps"))) + 
  mapview(add_links)

# P099270 CAREC 1B & 6B
segments1 <- c(1645, 1653, 1665, 1671, 1673, 1674, 1675, # Shymkent -> Khorgos (and beyond)
              1643, # Shymkent -> Tashkent
              1633, 1624, 1578, 1564, 1543) # Shymkent -> Aktobe

# P149952	East-West Highway Corridor Improvement (Georgia)
segments2 <- c(1364, 1337, 1328, 1317, 1309)

# P149955	Road Upgrading and Development Project (North Macedonia)
segments3 <- c(559, 560) # -> Segments are a bit longer than the actual projects!!

# P159220	Third/Fourth Phase of the Central Asia Regional Links Program (CARs-3) (Kyrgyz Republic)
#-> Seems not in the graph!!
segments4 <- c(1642, 1632, 1647)

# P173782	Kakheti Connectivity Improvement Project (Georgia)
segments5 <- c(1389, 1390, 1391) # 1391 goes outside the country to Baku

# P174379	Regional Connectivity and Development Project (Azerbaijan)
segments6 <- c(1454) # Way too long: project only goes from Salyan to Balsuvar. Maybe do 30% upgrade.

# P180153	Moldova Rural Connectivity Project (Moldowa)
segments7 <- c(776, 786, 791) # Not exactly project roads

# P500565	Transport Resilience and Connectivity Enhancement Project (Jezkazgan-Karagandy Section of TCITR (Middle Corridor))
segments8 <- c(1662) # Too long, only until Jezkazgan

segments <- c(segments1, segments2, segments3, segments4, segments5, segments6, segments7, segments8)

mapview(edges_real[segments, ])

# Saving 
graph_edges <- fread("data/transport_network/graph_orig.csv")
if(nrow(graph_edges) != nrow(edges_real)) stop("Mismatch!")

graph_edges$project <- FALSE
graph_edges$project[segments] <- TRUE

graph_edges |> fwrite("data/transport_network/graph_orig.csv")

# Plotting
pdf("figures/trans_ECA_network_projects.pdf", width = 10, height = 4.2)
tm_basemap("CartoDB.Positron", zoom = 5) +
  tm_shape(ss(edges_real, -segments)) +
  tm_lines(col = "grey70", lwd = 2) +
  tm_shape(ss(edges_real, segments)) +
  tm_lines(col = "orange", lwd = 2) +
  tm_shape(subset(nodes, !key_city)) + tm_dots(size = 0.1, fill = "grey50") +
  tm_shape(subset(nodes, key_city)) + tm_dots(size = 0.2) +
  tm_layout(frame = FALSE) 
dev.off()

pdf("figures/ECA_centroids_network_actual_discretized_middle_corridor_speed_EWTM.pdf", width = 17, height = 7)
tm_basemap("Esri.WorldTopoMap", zoom = 6) +
  tm_shape(edges_real_ext) +
  tm_lines(col = "speed_kmh",
           col.scale = tm_scale_continuous(5, values = "turbo"),
           col.legend = tm_legend("Speed in km/h", position = c("left", "top"), frame = FALSE, 
                                  text.size = 1, title.size = 1.2, margins = c(0, -0.5, 0, 0),
                                  title.padding = c(0, 0, -0.5, 0), 
                                  item.space = 0, item.height = 2, item.width = 0.5)) +
  tm_shape(subset(nodes, key_city)) + tm_dots(size = 0.2) +
  tm_shape(subset(nodes, !key_city)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE) #, inner.margins = c(0.1, 0.1, 0.1, 0.1))
dev.off()



# # Plot population and productivity (for GE Calibration) ------------------------------------------------------
# 
# graph_nodes <- fread("data/transport_network/csv/graph_nodes.csv") 
# graph_edges <- fread("data/transport_network/csv/graph_orig.csv") 
# 
# # Now: Plotting Productivity
# graph_nodes %<>%
#   mutate(citys = iif(product > 4L, 5L, product), 
#          prod_in = replace_na(37*outflows/population, 0),
#          prod = IWI + prod_in, 
#          total_prod = prod * population, 
#          total_dom_prod = IWI * population) %>%
#   st_as_sf(coords = c("lon", "lat"), crs = 4326)
# 
# # Check: This is imports divided by African GDP: 26.5%
# with(graph_nodes, sum(prod_in*population)/sum(IWI*population))
# with(graph_nodes, sum(prod*population)/sum(IWI*population))
# 
# table(graph_nodes$citys)
# 
# # <Table 8>
# graph_nodes |> qDT() |> 
#   collap(prod_in ~ citys, list(fsum, fmean, fmedian), give.names = FALSE) |> 
#   transpose(make.names = "citys", keep.names = "stat") |> 
#   tfmv(is.numeric, scales::label_number(scale_cut = scales::cut_short_scale(), accuracy = 0.01)) |>
#   xtable::xtable() |> print(include.r = FALSE, booktabs = TRUE)
# 
# # load("data/transport_network/trans_CEMAC_network.RData")
# edges_real <- qread("data/transport_network/edges_real_simplified.qs")
# 
# # <Figure 33>
# pdf("figures/trans_CEMAC_network_GE_parameterization_latest_MACR_90kmh_google.pdf", width = 9.5, height = 12)
# tm_basemap("CartoDB.Positron", zoom = 5) +
#   tm_shape(mutate(edges_real, speed_kmh = (edges$distance/1000)/(edges$duration/60^2)) |> 
#              rowbind(mutate(add_links, speed_kmh = 0), fill = TRUE)) +
#   tm_lines(col = "speed_kmh", 
#            col.legend = tm_legend("Speed (km/h)", position = tm_pos_in(0.57, 0.4), stack = "h", frame = FALSE, text.size = 1.2, title.size = 1.4),
#            col.scale = tm_scale_continuous(values = "turbo"), # 7, 
#            lwd = 2) + 
#   tm_shape(subset(graph_nodes, outflows > 0) |>
#              mutate(ofl = round(outflows / 1e6, 1))) + 
#   tm_dots(size = "ofl", 
#           size.scale = tm_scale_intervals(5, style = "jenks", values.scale = 2.5), # 
#           size.legend = tm_legend("Port Outflows (M)", position = tm_pos_in(0.57, 0.4), frame = FALSE, text.size = 1.2, title.size = 1.4),
#           fill = scales::alpha("black", 0.25)) +
#   tm_shape(subset(graph_nodes, population > 0) |> 
#              mutate(pop = population / 1000, prod = gdp / population)) + 
#   tm_dots(size = "pop", 
#           size.scale = tm_scale_intervals(breaks = c(0, 200, 1000, Inf),
#                                           values = c(1.5, 3, 5)*0.2),
#           size.legend = tm_legend("Population (K)", position = tm_pos_in(0.57, 0.17), stack = "h", frame = FALSE, text.size = 1.2, title.size = 1.4),
#           fill = "prod",
#           size.free = TRUE,
#           fill.scale = tm_scale_intervals(4, values = "inferno"), #viridis::inferno(5, alpha = 0.5, direction = -1)),
#           fill.legend = tm_legend("Productivity (GDP)", position = tm_pos_in(0.57, 0.17), frame = FALSE, text.size = 1.2, title.size = 1.4)) +
#   tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
#   tm_layout(frame = FALSE)
# dev.off()