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
nodes_param <- fread("data/transport_network/graph_nodes_MA_100kmh.csv")
nodes_param |> select(lon, lat) |> qM() |> subtract(st_coordinates(nodes)) |> descr() |> print(digits = 7)
nodes_param$product <- nodes$product
nodes_param$iso3c <- nodes$iso3c
# nodes_param |> atomic_elem() |> qDT() |> fwrite("data/transport_network/graph_nodes_MA_100kmh.csv")
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


# # Plot population and productivity (for GE Calibration) ------------------------------------------------------
# 
# graph_nodes <- fread("data/transport_network/csv/graph_nodes_MA_100kmh.csv") 
# graph_edges <- fread("data/transport_network/csv/graph_orig_MA_100kmh.csv") 
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