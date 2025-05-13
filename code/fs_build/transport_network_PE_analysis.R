#####################################################################
# ECA Transport Network: Partial Equilibrium Analysis
#####################################################################

library(fastverse)
set_collapse(mask = c("manip", "helper", "special"), nthreads = 4)
fastverse_extend(qs, sf, units, sfnetworks, tmap, install = TRUE)
source("code/fs_build/helpers.R")
fastverses()


# All essential objects from previous sections
edges_real <- qread("data/transport_network/edges_real_simplified.qs") # |> 
  # rmapshaper::ms_simplify(keep = 0.5)
if(!all.equal(edges$gravity_rd, edges_real$gravity_rd)) stop("Mismatch!")
tfm(edges_real) <- atomic_elem(edges)

# Average route efficiency per edge
aere <- unattrib(mean(edges$sp_distance / edges$distance)) 
print(aere)

# Adding distance to new edges based on average edge route efficiency
settfm(add_links, sp_distance = unattrib(st_length(geometry)))
settfm(add_links, distance = sp_distance / aere)


# Shortest Paths ----------------------------------------------------------

if(!net |> activate("edges") |> tidygraph::as_tibble() |> select(from, to) |> 
   atomic_elem() |> all.equal(atomic_elem(select(edges, from, to)))) stop("Mismatch") # Check

distances <- st_network_cost(net, weights = "distance") # distance in m
n_links <- igraph::distances(net)
range(n_links)

# Checks 
identical(st_geometry(net, "nodes"), nodes$geometry)
dist(unattrib(st_coordinates(nodes$geometry)), unattrib(qM(dist_ttime_mats$sources)))
sp_distances <- st_distance(nodes)

# Network Route Efficiency
nre <- mean(sp_distances / distances, na.rm = TRUE)
rnre <- mean(sp_distances / dist_ttime_mats$distances, na.rm = TRUE) # Real NRE

# Now adding edges
identical(st_geometry(net, "edges"), edges$geometry)
net_ext_data <- rbind(select(edges, distance, geometry), select(add_links, distance, geometry))
net_ext <- as_sfnetwork(net_ext_data, directed = FALSE)
plot(net_ext)
identical(st_geometry(net_ext, "nodes"), nodes$geometry) # Not the case, thus need to recalculate spherical distance as well
ind_ext <- ckmatch(nodes_coord, mctl(st_coordinates(st_geometry(net_ext, "nodes"))))
sp_distances_ext <- st_distance(st_geometry(net_ext, "nodes"))[ind_ext, ind_ext]
identical(sp_distances_ext, sp_distances)

distances_ext <- st_network_cost(net_ext, weights = "distance")[ind_ext, ind_ext]
all.equal(sp_distances_ext, sp_distances)
nre_ext <- mean(sp_distances_ext / distances_ext, na.rm = TRUE)

nre_ext / nre
rnre * (nre_ext / nre) # Reported increase
mean(distances / distances_ext, na.rm = TRUE) 

# Per link gain in NRE: takes a few mins
add_links$nre_per_link <- sapply(seq_row(add_links), function(i) {
  net_extd = as_sfnetwork(rbind(select(edges, distance), 
                                subset(add_links, i, distance)), directed = FALSE)
  distances_extd = st_network_cost(net_extd, weights = "distance")
  ind = ckmatch(nodes_coord, mctl(st_coordinates(st_geometry(net_extd, "nodes"))))
  mean(sp_distances / distances_extd[ind, ind], na.rm = TRUE)
})

# Percent increase
add_links$nre_gain_perc <- (unattrib(add_links$nre_per_link / nre) - 1) * 100
descr(add_links$nre_gain_perc)

# Plot percent increase
pdf("figures/PE/trans_ECA_network_NRE_gain_perc.pdf", width = 10, height = 4.2)
tm_basemap("CartoDB.Positron", zoom = 4) +
  tm_shape(edges_real) +
  tm_lines(lwd = 2, col = "grey70") +
  tm_shape(add_links) + 
  tm_lines(col = "nre_gain_perc", 
           col.scale = tm_scale_intervals(values = "turbo", breaks = c(0, 0.025, 0.1, 0.25, 0.5, Inf)),   
           col.legend = tm_legend(expression(Delta~"%"~"NRE"), position = c("right", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE) 
dev.off()


# Gravity weighted versions
gravity <- replace_inf(tcrossprod(nodes$population) / sp_distances, set = TRUE) |> unclass()
nre_wtd <- fmean(unattrib(sp_distances / distances), w = gravity)
rnre_wtd <- fmean(unattrib(sp_distances / dist_ttime_mats$distances), w = gravity)
all.equal(sp_distances_ext, sp_distances)
nre_ext_wtd <- fmean(unattrib(sp_distances_ext / distances_ext), w = gravity)

rnre_wtd
nre_ext_wtd / nre_wtd
rnre_wtd * (nre_ext_wtd / nre_wtd) # Reported increase

# Per link gain in NRE: weighted: takes a few mins
add_links$nre_wtd_per_link <- sapply(seq_row(add_links), function(i) {
  net_extd = as_sfnetwork(rbind(select(edges, distance), 
                                subset(add_links, i, distance)), directed = FALSE)
  distances_extd = st_network_cost(net_extd, weights = "distance")
  ind = ckmatch(nodes_coord, mctl(st_coordinates(st_geometry(net_extd, "nodes"))))
  fmean(unattrib(sp_distances / distances_extd[ind, ind]), w = gravity)
})
# Percent increase
add_links$nre_wtd_gain_perc <- (unattrib(add_links$nre_wtd_per_link / nre_wtd) - 1) * 100
descr(add_links$nre_wtd_gain_perc)

# Plot percent increase
# <Figure 19: RHS>
pdf("figures/PE/trans_ECA_network_NRE_wtd_gain_perc.pdf", width = 10, height = 4.2)
tm_basemap("CartoDB.Positron", zoom = 5) +
  tm_shape(edges_real) +
  tm_lines(lwd = 2, col = "grey70") +
  tm_shape(add_links) + 
  tm_lines(col = "nre_wtd_gain_perc", 
           col.scale = tm_scale_intervals(values = "turbo", breaks = c(0, 0.025, 0.1, 0.25, 0.5, Inf)),   
           col.legend = tm_legend(expression(Delta~"%"~"NRE WTD"), position = c("right", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE) 
dev.off()


# Market Access -------------------------------------------------------------------------------

# Matching to Locan GDP, and choosing distance-weighted nearest cell value within 50km
countries <- wbstats::wb_cachelist$countries %$% iso3c[region_iso3c %==% "ECS"] %>% 
  c("XKO") %>% setdiff(c("GRL", "ISL", "GBR", "IRL", "NOR", "SWE", "FIN", "FRO", "IMN", "CYP"))

# From geojson.io
ECA_outline <- st_read("/Users/sebastiankrantz/Documents/ECAOptimalTransport/data/outline/layers") |> 
  st_buffer(as_units(5, "km"))

# Local GDP data: 0.5 degree, censored for density below 0.01
LGDP <- fread("/Users/sebastiankrantz/Documents/Data/LocalGDP/0_25deg/final_GDP_0_25deg_postadjust_pop_dens_0_01_adjust.csv") |> 
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
gvr(LGDP, "GDPC") %*=% 1e6

pdf("figures/PE/cell_GDPC_const_2017_PPP_ECA_GRID.pdf", width = 10, height = 4.2)
tm_basemap("CartoDB.Positron", zoom = 5) +
  tm_shape(edges_real) + tm_lines(col = "grey30") + 
  tm_shape(st_as_sf(LGDP, coords = c("longitude", "latitude"), crs = 4326)) + 
  tm_symbols(fill = "cell_GDPC_const_2017_PPP", size = 0.235, shape = 15, fill_alpha = 0.6, # ?points
             fill.scale = tm_scale_intervals(7, style = "fisher", values = "inferno"),
             fill.legend = tm_legend("Cell GDPC", position = c("left", "top"), frame = FALSE, 
                                     text.size = 1, title.size = 1.2, 
                                     title.padding = c(-0.3, 0, -0.5, 0), 
                                     item.space = 0)) + 
  tm_layout(frame = FALSE)
dev.off()

outcomes_coords <- LGDP |> select(longitude, latitude) |> qM()
nodes_coord_mat <- st_coordinates(nodes) # graph_nodes |> select(lon, lat) |> qM() 
nodes$N_GDPC <- nodes$GDPC <- NA_real_
for (i in seq_row(nodes)) {
  d = geodist::geodist(nodes_coord_mat[i, , drop = FALSE], outcomes_coords, measure = "haversine") 
  ind = which(d < 50e3)
  if(!length(ind)) ind = which(d < min(d) + 50e3)
  w = 1 / d[ind]
  nodes$GDPC[i] = fmean.default(LGDP$cell_GDPC_const_2017_PPP[ind], w = w)
  nodes$N_GDPC[i] = length(ind)
}
settfm(nodes, GDP = GDPC * population)
rm(d, ind, w, nodes_coord_mat, outcomes_coords, GDPC); gc()
fndistinct(atomic_elem(nodes))

pdf("figures/PE/trans_ECA_network_cell_GDPC_const_2017_PPP.pdf", width = 10, height = 4.2)
tm_basemap("CartoDB.Positron", zoom = 5) +
  tm_shape(edges_real) + tm_lines(col = "grey30") + 
  tm_shape(nodes) + 
  tm_dots(fill = "GDPC", size = 0.25,
          fill.scale = tm_scale_intervals(7, style = "pretty", values = "inferno"),
          fill.legend = tm_legend("Productivity (GDPC)", position = c("left", "top"), frame = FALSE, 
                                  text.size = 1, title.size = 1.2, margins = c(0, -0.5, 0, 0),
                                  title.padding = c(0, 0, -0.5, 0), 
                                  item.space = 0)) + 
  tm_layout(frame = FALSE)
dev.off()


# Computing total market access ----------------------------------------------------

(MA_real <- total_MA(dist_ttime_mats$distances, nodes$GDP))
(MA <- total_MA(distances, nodes$GDP)) # distances^3.8

# Total gain
(MA_ext <- total_MA(distances_ext, nodes$GDP)) # distances^3.8 

MA_ext / MA

# Needed for later
ma_gain_per_km <- (MA_ext - MA) * 1000

# Compute change in MA from each link
add_links$MA_per_link <- sapply(seq_row(add_links), function(i) {
  nete = as_sfnetwork(rbind(select(edges, distance), 
                            subset(add_links, i, distance)), directed = FALSE)
  inv_dist = 1 / unclass(st_network_cost(nete, weights = "distance"))
  diag(inv_dist) = 0
  ind = ckmatch(mctl(st_coordinates(st_geometry(nete, "nodes"))), nodes_coord)
  sum(inv_dist %*% nodes$GDP[ind])
})
# Percent increase
add_links$MA_gain_perc <- (add_links$MA_per_link / MA - 1) * 100
descr(add_links$MA_gain_perc)

# <Figure 20: LHS> (Use distances^3.8 above to generate RHS)
pdf("figures/PE/trans_ECA_network_MA_gain_perc.pdf", width = 10, height = 4.2)
tm_basemap("CartoDB.Positron", zoom = 5) +
  tm_shape(edges_real) +
  tm_lines(lwd = 2, col = "grey70") +
  tm_shape(add_links) + 
  tm_lines(col = "MA_gain_perc", 
           col.scale = tm_scale_intervals(values = "turbo", breaks = c(0, 0.025, 0.1, 0.25, 0.5, Inf)),
           col.legend = tm_legend(expression(Delta~"%"~"MA"), position = c("left", "top"), frame = FALSE, 
                                  text.size = 1, title.size = 1.2, margins = c(0, -0.5, 0, 0),
                                  title.padding = c(0, 0, -0.5, 0), 
                                  item.space = 0), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE) 
dev.off()



# Estimating Network Building Costs -------------------------------------------

# This commented code shows how the 3000m buffer around links is computed and applied
add_links_buff_3km <- st_buffer(add_links, as_units(3000, "m"))
edges_buff_3km <- st_buffer(edges_real, as_units(3000, "m"))

# Adding Ruggedness: https://diegopuga.org/data/rugged/
rugg <- terra::rast("/Users/sebastiankrantz/Documents/Data/Ruggedness/tri.txt")
rugg_area <- terra::rast("/Users/sebastiankrantz/Documents/Data/Ruggedness/cellarea.txt")
add_links$rugg <- exactextractr::exact_extract(rugg, add_links_buff_3km, weights = rugg_area, fun = "weighted_mean")
edges$rugg <- exactextractr::exact_extract(rugg, edges_buff_3km, weights = rugg_area, fun = "weighted_mean")
rm(rugg, rugg_area)

# Adding Population (WorldPop 2020 1km2 global)
pop_wpop <- terra::rast("/Users/sebastiankrantz/Documents/Data/WorldPop/ppp_2020_1km_Aggregated.tif")
add_links$pop_wpop <- exactextractr::exact_extract(pop_wpop, add_links_buff_3km, fun = "sum")
add_links$pop_wpop_km2 <- unattrib(add_links$pop_wpop / (st_area(add_links_buff_3km) / 1e6))
edges$pop_wpop <- exactextractr::exact_extract(pop_wpop, edges_buff_3km, fun = "sum")
edges$pop_wpop_km2 <- unattrib(edges$pop_wpop / (st_area(edges_buff_3km) / 1e6))
# Cleanup
rm(add_links_buff_3km, edges_buff_3km, pop_wpop); gc()
all.equal(select(qDF(edges_real), from, to), select(qDF(edges), from, to))
tfm(edges_real) <- atomic_elem(edges)

# Plot Ruggedness
# <Figure 23: LHS>
pdf("figures/PE/trans_ECA_network_rugg.pdf", width = 10, height = 4.2)
tm_basemap("Esri.WorldTopoMap", zoom = 5) +
  tm_shape(mutate(rbind(select(edges_real, rugg), select(add_links, rugg)), rugg = rugg / 1000)) +
  tm_lines(col = "rugg",
           col.scale = tm_scale_continuous_log1p(7, values = "turbo"),
           col.legend = tm_legend("Ruggedness", position = c("left", "top"), frame = FALSE, 
                                  text.size = 1, title.size = 1.2, margins = c(0, -0.5, 0, 0),
                                  title.padding = c(0, 0, -0.5, 0), 
                                  item.space = 0, item.height = 1.5, item.width = 0.5), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE)
dev.off()

add_links |> gvr("_km2") |> descr()
edges |> gvr("_km2") |> descr()

# Plot Population Density
# <Figure 23: RHS>
pdf("figures/PE/trans_ECA_network_pop_wpop_km2.pdf", width = 10, height = 4.2)
tm_basemap("CartoDB.Positron", zoom = 5) +
  tm_shape(rbind(select(edges_real, pop_wpop_km2), select(add_links, pop_wpop_km2))) +
  tm_lines(col = "pop_wpop_km2",
           col.scale = tm_scale_continuous_log1p(7, values = "turbo"),
           col.legend = tm_legend(expression("Population/km"^2), position = c("left", "top"), frame = FALSE, 
                                  text.size = 1, title.size = 1.2, margins = c(0, -0.5, 0, 0),
                                  title.padding = c(0, 0, -0.5, 0), 
                                  item.space = 0, item.height = 1.5, item.width = 0.5), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE)
dev.off()



# Estimating Road Costs following my Optimal Roads Paper: ------------------------------

# Estimates from Collier, Kirchberger & Söderbom (2016)
# Pop Coef
mean(c(0.11, 0.088, 0.082,   # Table 4
       0.077, 0.083, 0.074)) # Taable 5

# Calibrating cost eauation to match median 2L Highway construction cost in Africa (611 million/km)
mean(with(edges, exp(log(100e3) - 0.11 * (distance > 50e3) + 0.12 * log(rugg) + 0.085 * log(pop_wpop_km2+1)))) / 1000 

# Applying Equation
settfm(edges, cost_km = exp(log(100e3) - 0.11 * (distance > 50e3) + 0.12 * log(rugg) + 0.085 * log(pop_wpop_km2+1))) # 0.00085 * pop_wpop_km2
descr(edges$cost_km)

settfm(add_links, cost_km = exp(log(100e3) - 0.11 * (distance > 50e3) + 0.12 * log(rugg) + 0.085 * log(pop_wpop_km2+1))) # 0.00085 * pop_wpop_km2
descr(add_links$cost_km)

descr(rbind(select(edges, cost_km), select(add_links, cost_km)))

# Total Network Length and Cost
sum(edges$distance / 1000) / 1e3
sum(edges$cost_km * edges$distance / 1000) / 1e9
sum(add_links$distance / 1000) / 1e3
sum(add_links$cost_km * add_links$distance / 1000) / 1e9

edges_real$cost_km = edges$cost_km

# Plots
# <Figure A12: LHS>
pdf("figures/PE/trans_ECA_network_cost_km.pdf", width = 10, height = 4.2)
tm_basemap("CartoDB.Positron", zoom = 5) +
  tm_shape(rbind(select(edges_real, cost_km), select(add_links, cost_km))) +
  tm_lines(col = "cost_km",
           col.scale = tm_scale_continuous(7, values = "turbo"), 
           col.legend = tm_legend("USD'15/km", position = c("left", "top"), frame = FALSE, 
                                  text.size = 1, title.size = 1.2, margins = c(0, -0.5, 0, 0),
                                  title.padding = c(0, 0, -0.5, 0), 
                                  item.space = 0, item.height = 1.5, item.width = 0.5), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE)
dev.off()


# Cost-Benefit Analysis: Adding New Links -------------------------------------------

# Total Cost-Benefit Ratios

ma_gain_per_km / 1e9 # MA gain in billions
ma_gain_per_km / sum(with(add_links, cost_km * distance / 1000)) # MA gain per investment

# MA Gain per Dollar
settfm(add_links, MA_gain_pusd = perch_to_diff(MA_per_link, MA_gain_perc) * 1000 / (cost_km * distance / 1000)) # * 1216
descr(add_links$MA_gain_pusd)
proportions(table(add_links$MA_gain_pusd < 0.5))

# <Figure 24: LHS (Top)>
pdf("figures/PE/trans_ECA_network_MA_gain_pusd.pdf", width = 10, height = 4.2)
tm_basemap("CartoDB.Positron", zoom = 5) +
  tm_shape(edges_real) +
  tm_lines(lwd = 2, col = "grey70") +
  tm_shape(add_links) + 
  tm_lines(col = "MA_gain_pusd", 
           col.scale = tm_scale_intervals(values = "turbo", breaks = c(0, 0.1, 0.2, 0.5, 1, 2, Inf)/10),
           col.legend = tm_legend(expression(Delta~"MA/USD"), position = c("left", "top"), frame = FALSE, 
                                  text.size = 1, title.size = 1.2, margins = c(0, -0.5, 0, 0),
                                  title.padding = c(0, 0, -0.5, 0), 
                                  item.space = 0, item.height = 1.5, item.width = 0.5), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE) 
dev.off()

# Consensus Package
settfm(add_links, consensus = MA_gain_pusd > 0.01) 

# <Figure 24: RHS (Bottom)>
pdf("figures/PE/trans_ECA_network_MA_gain_pusd_cons.pdf", width = 10, height = 4.2)
tm_basemap("CartoDB.Positron", zoom = 5) +
  tm_shape(edges_real) +
  tm_lines(lwd = 2, col = "grey70") +
  tm_shape(subset(add_links, consensus, MA_gain_pusd)) + 
  tm_lines(col = "MA_gain_pusd", 
           col.scale = tm_scale_intervals(values = "turbo", breaks = c(0, 0.1, 0.2, 0.5, 1, 2, Inf)/10),
           col.legend = tm_legend(expression(Delta~"MA/USD"), position = c("left", "top"), frame = FALSE, 
                                  text.size = 1, title.size = 1.2, margins = c(0, -0.5, 0, 0),
                                  title.padding = c(0, 0, -0.5, 0), 
                                  item.space = 0, item.height = 1.5, item.width = 0.5), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE) 
dev.off()

# Consensus Gains
nrow(subset(add_links, consensus)) / nrow(add_links)
subset(add_links, consensus) |> with(sum(cost_km * distance / 1000)) |> divide_by(1e9)

net_ext_cons <- as_sfnetwork(rbind(select(edges, distance), 
                                   subset(add_links, consensus, distance)), directed = FALSE)
plot(net_ext_cons)
ind_ext_cons <- ckmatch(nodes_coord, mctl(st_coordinates(st_geometry(net_ext_cons, "nodes"))))
identical(st_distance(st_geometry(net_ext_cons, "nodes"))[ind_ext_cons, ind_ext_cons], sp_distances)

distances_ext_cons <- st_network_cost(net_ext_cons, weights = "distance")[ind_ext_cons, ind_ext_cons]

# Total gain
(MA_ext_cons <- total_MA(distances_ext_cons, nodes$GDP)) 

MA_ext_cons / MA

ma_gain_per_km_cons <- (MA_ext_cons - MA) * 1000
ma_gain_per_km_cons / 1e9 # MA gain in billions
ma_gain_per_km_cons / sum(with(subset(add_links, consensus), cost_km * distance / 1000)) # MA gain per investment


# Now: Improving Existing Links ------------------------------------------------------------------------

settfm(edges, speed_kmh = (distance / 1000) / (duration / 60))
descr(edges$speed_kmh)

# <Figure 25: LHS>
hist(edges$speed_kmh, breaks = 80, xlab = "Average Link Speed in km/h", main = NULL)
dev.copy(pdf, "figures/PE/trans_ECA_network_average_link_speed_hist.pdf", width = 8, height = 7)
dev.off()

# Inspect
# edges |> select(speed_kmh) |> mapview::mapview() 

tfm(edges_real) <- atomic_elem(edges)

# Plot
# <Figure 25: RHS>
pdf("figures/PE/trans_ECA_network_average_link_speed.pdf", width = 10, height = 4.2)
tm_basemap("CartoDB.Positron", zoom = 5) +
  tm_shape(edges_real) +
  tm_lines(col = "speed_kmh", 
           col.scale = tm_scale_continuous(5, values = "turbo"),
           col.legend = tm_legend("Speed in km/h", position = c("left", "top"), frame = FALSE, 
                                  text.size = 1, title.size = 1.2, margins = c(0, -0.5, 0, 0),
                                  title.padding = c(0, 0, -0.5, 0), 
                                  item.space = 0, item.height = 2.5, item.width = 0.5), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE) 
dev.off()

# Computing Times 
if(!net |> activate("edges") |> tidygraph::as_tibble() |> select(from, to) |> 
   atomic_elem() |> all.equal(atomic_elem(select(edges, from, to)))) stop("Mismatch") # Check

times <- st_network_cost(net, weights = edges$duration)
edges %<>% mutate(speed_kmh_imp = iif(speed_kmh < 100, 100, speed_kmh),
                  duration_imp = duration * speed_kmh / speed_kmh_imp)
descr(edges, cols = .c(duration, duration_imp))
times_imp <- st_network_cost(net, weights = edges$duration_imp)

# Computing total real market access
(MA_real <- total_MA(dist_ttime_mats$durations, nodes$GDP)) # Original: 1748.128 billion USD/min
(MA <- total_MA(times, nodes$GDP)) #  

# Total gain
(MA_imp <- total_MA(times_imp, nodes$GDP))

MA_imp / MA
# Gain from original: 
(MA_imp / MA) * MA_real

# Needed for later
ma_gain_per_min <- MA_imp - MA

edges$MA_100_min_speed <- sapply(seq_row(edges), function(i) {
  w = copyv(edges$duration, i, edges$duration_imp, vind1 = TRUE)
  inv_dur = 1 / unclass(st_network_cost(net, weights = w))
  diag(inv_dur) = 0 
  sum(inv_dur %*% nodes$GDP) 
})
# Percent increase
edges$MA_100_min_speed_perc <- (edges$MA_100_min_speed / MA - 1) * 100
descr(edges$MA_100_min_speed_perc)

tfm(edges_real) <- atomic_elem(edges)

# <Figure 26: A>
pdf("figures/PE/trans_ECA_network_MA_100_min_speed_perc.pdf", width = 10, height = 4.2)
tm_basemap("CartoDB.Positron", zoom = 5) +
  tm_shape(edges_real) +
  tm_lines(col = "MA_100_min_speed_perc", 
           col.scale = tm_scale_intervals(values = "brewer.yl_or_rd", breaks = c(0, 0.01, 0.025, 0.1, 0.25, 0.5, Inf)/5),
           col.legend = tm_legend(expression(Delta~"%"~"MA [GDP/min]"),
                                  position = c("left", "top"), frame = FALSE, 
                                  text.size = 1, title.size = 1.2, margins = c(0, -0.5, 0, 0),
                                  title.padding = c(0, 0, -0.5, 0), 
                                  item.space = 0, item.height = 1.5, item.width = 0.5), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE)
dev.off()

# Considering the addition of proposed links under 90km/h or 65km/h assumption
settfm(add_links, 
       duration_100kmh = distance / kmh_to_mmin(90), 
       duration_65kmh = distance / kmh_to_mmin(65))

# Temporary networks as needed
net_ext_tmp <- as_sfnetwork(rbind(select(edges, duration = duration_imp), 
                                  select(add_links, duration = duration_100kmh)), directed = FALSE)
ind_ext_tmp <- ckmatch(nodes_coord, mctl(st_coordinates(st_geometry(net_ext_tmp, "nodes"))))
times_ext_tmp <- st_network_cost(net_ext_tmp, weights = "duration")[ind_ext_tmp, ind_ext_tmp]

# Total gain
(MA_tmp <- total_MA(times_ext_tmp, nodes$GDP)) / 1e9

MA_tmp / MA # * MA_real
(MA_tmp - MA) / 1e9
rm(list = ls()[endsWith(ls(), "_tmp")]); gc()


# Link Upgrading Costs ---------------------------------------------------

settfm(edges, upgrade_cat = nif(speed_kmh < 60, "Upgrade: <60km/h", speed_kmh >= 60 & speed_kmh < 80, "Mixed Works: 60-80km/h", 
                                speed_kmh >= 80 & speed_kmh < 100, "Resurfacing: 80-100km/h", speed_kmh >= 100, "Nothing: >100km/h") |> 
         factor(levels = c("Nothing: >100km/h", "Resurfacing: 80-100km/h", "Mixed Works: 60-80km/h", "Upgrade: <60km/h")))
table(edges$upgrade_cat)
anyNA(edges$upgrade_cat)

# <Figure 29: LHS>
tfm(edges_real) <- atomic_elem(edges)
pdf("figures/PE/trans_ECA_network_type_of_work.pdf", width = 10, height = 4.2)
tm_basemap("CartoDB.Positron", zoom = 5) +
  tm_shape(edges_real) +
  tm_lines(col = "upgrade_cat", 
           col.scale = tm_scale_categorical(values = "turbo"),
           col.legend = tm_legend("Type of Work", position = c("left", "top"), frame = FALSE, 
                                  text.size = 1, title.size = 1.2, margins = c(0, -0.5, 0, 0),
                                  title.padding = c(0, 0, -0.5, 0), 
                                  item.space = 0, item.height = 1.5, item.width = 0.5), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE)
dev.off()

# Now costing the categories
descr(with(edges, # subset(edges, upgrade_cat == "Resurfacing"), 
           exp(log(28.4e3 * 100/120) - 0.11 * (distance > 50e3) + 0.12 * log(rugg) + 0.085 * log(pop_wpop_km2+1)))) 
descr(with(edges, # subset(edges, upgrade_cat == "Mixed Works"), 
           exp(log(64.6e3 * 100/120) - 0.11 * (distance > 50e3) + 0.12 * log(rugg) + 0.085 * log(pop_wpop_km2+1)))) 
descr(with(edges, # subset(edges, upgrade_cat == "Upgrade"), 
           exp(log(101.6e3 * 100/120) - 0.11 * (distance > 50e3) + 0.12 * log(rugg) + 0.085 * log(pop_wpop_km2+1)))) 

edges %<>%
  mutate(ug_cost_km = -0.11 * (distance > 50e3) + 0.12 * log(rugg) + 0.085 * log(pop_wpop_km2+1), 
         ug_cost_km = nif(upgrade_cat == "Resurfacing: 80-100km/h", clip5perc(exp(ug_cost_km + log(28.4e3 * 100/120))),
                          upgrade_cat == "Mixed Works: 60-80km/h", clip5perc(exp(ug_cost_km + log(64.6e3 * 100/120))),
                          upgrade_cat == "Upgrade: <60km/h", clip5perc(exp(ug_cost_km + log(101.6e3 * 100/120))), 
                          upgrade_cat == "Nothing: >100km/h", 0)) 

descr(edges$ug_cost_km)
descr(edges, ug_cost_km ~ upgrade_cat)
hist(edges$ug_cost_km, breaks = 80)

# <Figure 29: RHS>
tfm(edges_real) <- atomic_elem(edges)
pdf("figures/PE/trans_ECA_network_upgrading_costs.pdf", width = 10, height = 4.2)
tm_basemap("CartoDB.Positron", zoom = 5) +
  tm_shape(edges_real) +
  tm_lines(col = "ug_cost_km", 
           col.scale = tm_scale_continuous(7, values = "brewer.yl_or_rd"), 
           col.legend = tm_legend("USD'15/km", position = c("left", "top"), frame = FALSE, 
                                  text.size = 1, title.size = 1.2, margins = c(0, -0.5, 0, 0),
                                  title.padding = c(0, 0, -0.5, 0), 
                                  item.space = 0, item.height = 1.5, item.width = 0.5), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE)
dev.off()

# Total Costs and Breakdown
sum(edges$ug_cost_km * edges$distance / 1000) / 1e9
fsum(edges$ug_cost_km * edges$distance / 1000, edges$upgrade_cat) / 1e9


# Cost-Benefit Analysis ---------------------------------------------------

ma_gain_per_min / 1e9 # MA gain in billions
ma_gain_per_min / sum(with(edges, ug_cost_km * distance / 1000)) # MA gain per investment

# MA Gain per Dollar
settfm(edges, MA_gain_pusd = perch_to_diff(MA_100_min_speed, MA_100_min_speed_perc) / (ug_cost_km * distance / 1000)) # * 1216
edges$MA_gain_pusd |> replace_inf(set = TRUE)
# edges$MA_gain_pusd |> replace_na(0, set = TRUE)
descr(edges$MA_gain_pusd)
proportions(table(edges$MA_gain_pusd < 0.5))
proportions(table(edges$MA_gain_pusd < 0.25))

# <Figure 30: LHS (Top)>
tfm(edges_real) <- atomic_elem(edges)
pdf("figures/PE/trans_ECA_network_MA_gain_100_min_speed_pusd.pdf", width = 10, height = 4.2)
tm_basemap("CartoDB.Positron", zoom = 5) +
  tm_shape(edges_real) + 
  tm_lines(col = "MA_gain_pusd", 
           col.scale = tm_scale_intervals(values = "-inferno", breaks = c(0, 0.1, 0.2, 0.5, 1, 2, 5, Inf)/10),
           col.legend = tm_legend(expression(Delta~"MA/USD"), position = c("left", "top"), frame = FALSE, 
                                  text.size = 1, title.size = 1.2, margins = c(0, -0.5, 0, 0),
                                  title.padding = c(0, 0, -0.5, 0), 
                                  item.space = 0, item.height = 1.5, item.width = 0.5), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE)
dev.off()

# Consensus Package
settfm(edges, consensus = is.finite(MA_gain_pusd) & MA_gain_pusd > 0.1)

# <Figure 30: RHS (Bottom)>
tfm(edges_real) <- atomic_elem(edges)
pdf("figures/PE/trans_ECA_network_MA_gain_100_min_speed_pusd_cons.pdf", width = 10, height = 4.2)
tm_basemap("CartoDB.Positron", zoom = 5) +
  tm_shape(subset(edges_real, !consensus)) + tm_lines(lwd = 2, col = "grey70") +
  tm_shape(subset(edges_real, consensus, MA_gain_pusd)) +
  tm_lines(col = "MA_gain_pusd", 
           col.scale = tm_scale_intervals(values = "turbo", breaks = c(0, 0.1, 0.2, 0.5, 1, 2, 5, Inf)/10),
           col.legend = tm_legend(expression(Delta~"MA/USD"), position = c("left", "top"), frame = FALSE, 
                                  text.size = 1, title.size = 1.2, margins = c(0, -0.5, 0, 0),
                                  title.padding = c(0, 0, -0.5, 0), 
                                  item.space = 0, item.height = 1.5, item.width = 0.5), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE)
dev.off()

# Consensus Gains
nrow(subset(edges, consensus)) / nrow(edges)
subset(edges, consensus) |> with(sum(ug_cost_km * distance / 1000)) |> divide_by(1e9)

net_imp_cons <- as_sfnetwork(rbind(subset(edges, !consensus, duration), 
                                   subset(edges, consensus, duration = duration_imp)), 
                             directed = FALSE)
# plot(net_imp_cons)
ind_imp_cons <- ckmatch(nodes_coord, mctl(st_coordinates(st_geometry(net_imp_cons, "nodes"))))
identical(st_distance(st_geometry(net_imp_cons, "nodes"))[ind_imp_cons, ind_imp_cons], sp_distances)

times_imp_cons <- st_network_cost(net_imp_cons, weights = "duration")[ind_imp_cons, ind_imp_cons]

# Total gain
MA_imp_cons <- total_MA(times_imp_cons, nodes$GDP)
MA_imp_cons / MA

ma_gain_per_min_cons <- MA_imp_cons - MA

ma_gain_per_min_cons / 1e9 # MA gain in billions
ma_gain_per_min_cons / sum(with(subset(edges, consensus), ug_cost_km * distance / 1000)) # MA gain per investment


# Cost-Benefit Analysis: Joint Scenarios ---------------------------------------------------

# Plot All Costs
all_costs <- rbind(existing = select(edges_real, cost_km = ug_cost_km, distance, duration, duration_imp), 
                     new = select(add_links, cost_km = cost_km, distance, duration_imp = duration_100kmh) |> 
                           transform(duration = duration_imp))
all_costs$type <- iif(rownames(all_costs) %ilike% "existing", "existing", "new")
descr(all_costs$cost_km)

# <Figure 31: LHS>
hist(all_costs$cost_km / 1000, breaks = 80, xlab = "Cost per Km in Thousands of 2015 USD", main = NULL)
dev.copy(pdf, "figures/PE/trans_ECA_network_all_costs_hist.pdf", width = 8, height = 7)
dev.off()

# <Figure 31: RHS>
pdf("figures/PE/trans_ECA_network_all_costs.pdf", width = 10, height = 4.2)
tm_basemap("CartoDB.Positron", zoom = 5) +
  tm_shape(all_costs) +
  tm_lines(col = "cost_km", 
           col.scale = tm_scale_continuous(7, values = "turbo"), 
           col.legend = tm_legend("USD'15/km", position = c("left", "top"), frame = FALSE, 
                                  text.size = 1, title.size = 1.2, margins = c(0, -0.5, 0, 0),
                                  title.padding = c(0, 0, -0.5, 0), 
                                  item.space = 0, item.height = 1.5, item.width = 0.5), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE)
dev.off()

# Total Costs
sum(all_costs$cost_km * all_costs$distance / 1000) / 1e9

# Need to repeat simulations for new links with MA denominated in time 

# Baseline Scenario
add_links$MA_per_link_100kmh <- sapply(seq_row(add_links), function(i) {
  nete = as_sfnetwork(rbind(select(edges, duration), 
                            subset(add_links, i, duration = duration_100kmh)), directed = FALSE)
  ind = ckmatch(mctl(st_coordinates(st_geometry(nete, "nodes"))), nodes_coord)
  inv_dur = 1 / unclass(st_network_cost(nete, weights = "duration"))
  diag(inv_dur) = 0
  sum(inv_dur %*% nodes$GDP[ind])
})
add_links$MA_per_link_100kmh_perc <- (add_links$MA_per_link_100kmh / MA - 1) * 100
descr(add_links$MA_per_link_100kmh_perc)

# Computing Cost-Benefit Ratios
settfm(add_links, 
   MA_gain_100kmh_pusd = perch_to_diff(MA_per_link_100kmh, MA_per_link_100kmh_perc) / (cost_km * distance / 1000) |> replace_inf()
)

# Combining Datasets
all_cb_ratios <- rbind(edges_real |> select(MA_gain_pusd), # MA_gain_pusd_bt_opt
                       add_links |> select(MA_gain_100kmh_pusd) |> rm_stub("100kmh_", regex = TRUE)) # MA_gain_100kmh_pusd_bt_opt
tfm(all_cb_ratios) <- all_costs |> atomic_elem()
descr(all_cb_ratios)

# <Figure 32: First 3 Plots>
pdf("figures/PE/trans_ECA_network_MA_gain_all_100kmh_pusd_inferno.pdf", width = 10, height = 4.2)
  tm_basemap("CartoDB.Positron", zoom = 5) +
      tm_shape(all_cb_ratios) + 
      tm_lines(col = "MA_gain_pusd", 
               col.scale = tm_scale_intervals(values = "-inferno", breaks = c(0, 0.1, 0.2, 0.5, 1, 2, 5, Inf)/5),
               col.legend = tm_legend(expression(Delta~"MA/USD"),
                                      position = c("left", "top"), frame = FALSE, 
                                      text.size = 1, title.size = 1.2, margins = c(0, -0.5, 0, 0),
                                      title.padding = c(0, 0, -0.5, 0), 
                                      item.space = 0, item.height = 1.5, item.width = 0.5), lwd = 2) + 
      tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
      tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
      tm_layout(frame = FALSE)
dev.off()



# Need again without real edges (for total gain evaluation)
all_cb_ratios_se <- rbind(edges |> select(MA_gain_pusd), 
                          add_links |> select(MA_gain_100kmh_pusd) |> 
                            rm_stub("100kmh_", regex = TRUE)) 
tfm(all_cb_ratios_se) <- all_costs |> atomic_elem()
descr(all_cb_ratios_se)


for (i in c(0.1, 0.25, 0.5)) {
  cat("MA Gain greatern than ", i, fill = TRUE)
  
  settfm(all_cb_ratios, 
         consensus = is.finite(MA_gain_pusd) & MA_gain_pusd > i, 
         MA_gain_pusd_cons = MA_gain_pusd)
  tfm(all_cb_ratios_se) <- atomic_elem(all_cb_ratios)
  
  # <Figure 32: Last 3 Plots>
  pdf(sprintf("figures/PE/trans_ECA_network_MA_gain_all_100kmh_pusd_cons_MAg%g.pdf", i), width = 10, height = 4.2)
  print(tm_basemap("CartoDB.Positron", zoom = 5) +
    tm_shape(subset(all_cb_ratios, !consensus)) + tm_lines(lwd = 2, col = "grey70") +
    tm_shape(subset(all_cb_ratios, consensus, MA_gain_pusd_cons)) +
    tm_lines(col = "MA_gain_pusd_cons", 
             col.scale = tm_scale_intervals(values = "turbo", breaks = c(0, 0.1, 0.2, 0.5, 1, 2, 5, Inf)/5),
             col.legend = tm_legend(expression(Delta~"MA/USD"), position = c("left", "top"), frame = FALSE, 
                                    text.size = 1, title.size = 1.2, margins = c(0, -0.5, 0, 0),
                                    title.padding = c(0, 0, -0.5, 0), 
                                    item.space = 0, item.height = 1.5, item.width = 0.5), lwd = 2) + 
    tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
    tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
    tm_layout(frame = FALSE))
  dev.off()
  
  Sys.sleep(2)
  
  # Consensus Gains
  all_cb_ratios %$% table(type, consensus) |> proportions(1)
  all_cb_ratios %$% table(type, consensus, w = distance) |> # addmargins() |> 
    t() |> fsum.matrix(TRA = "%")
  sum(subset(all_cb_ratios, consensus)$distance) / sum(all_cb_ratios$distance)
  
  # Cost
  subset(all_cb_ratios, consensus) |> with(sum(cost_km * distance / 1000)) |> divide_by(1e9) |> print()
  
  net_imp_cons <- as_sfnetwork(rbind(
    subset(all_cb_ratios_se, consensus, duration = duration_imp),
    subset(all_cb_ratios_se, !consensus & type == "existing", duration)), directed = FALSE)
  
  # plot(net_imp_cons)
  ind_imp_cons <- ckmatch(nodes_coord, mctl(st_coordinates(st_geometry(net_imp_cons, "nodes"))))
  identical(st_distance(st_geometry(net_imp_cons, "nodes"))[ind_imp_cons, ind_imp_cons], sp_distances)
  times_imp_cons <- st_network_cost(net_imp_cons, weights = "duration")[ind_imp_cons, ind_imp_cons]

  # Total gain: No Frictions
  MA_imp_cons <- total_MA(times_imp_cons, nodes$GDP) 
  print(MA_imp_cons / MA) 
  ma_gain_per_min_cons <- MA_imp_cons - MA 
  print(ma_gain_per_min_cons / 1e9) # MA gain in billions
  # MA gain per investment:
  print(ma_gain_per_min_cons / sum(with(subset(all_cb_ratios, consensus), cost_km * distance / 1000)))
}





#############################################
# Saving PE Results
#############################################


list(nodes = nodes,
     edges = edges, 
     add_links = add_links) |>
  qsave("results/transport_network/PE/PE_results_MA_100kmh.qs")

nodes |> transform(set_names(mctl(st_coordinates(geometry)), c("lon", "lat"))) |> 
  atomic_elem() |> qDT() |> fwrite("results/transport_network/PE/PE_results_nodes_MA_100kmh.csv")
edges |> atomic_elem() |> qDT() |> fwrite("results/transport_network/PE/PE_results_edges_MA_100kmh.csv")
add_links |> atomic_elem() |> qDT() |> fwrite("results/transport_network/PE/PE_results_add_links_MA_100kmh.csv")

# Also saving Simplified real edges
edges_real <- qread("data/transport_network/edges_real_simplified.qs")
all.equal(edges_real$gravity_rd, edges$gravity_rd)
tfm(edges_real) <- atomic_elem(edges)

edges_real |> qsave("results/transport_network/PE/PE_edges_real_MA_100kmh.qs")


#############################################
# Saving Graphs for OTN (GE Analysis)
#############################################

qsu(edges)

graph_orig <- edges |> qDT() |> 
  select(from, to, sp_distance, distance, duration, speed_kmh, 
         speed_kmh_imp, duration_imp, rugg, pop_wpop, pop_wpop_km2, cost_km, upgrade_cat, ug_cost_km)

graph_add <- add_links |> qDT() |> 
  select(from, to, sp_distance, distance, duration_100kmh, duration_65kmh, rugg, pop_wpop, pop_wpop_km2, cost_km)

graph_nodes <- nodes |> qDT() |> 
  transform(set_names(mctl(st_coordinates(geometry)), c("lon", "lat"))) |> 
  select(lon, lat, city_ascii, key_city, population, GDPC, N_GDPC, GDP) 


# # Consistency Checks: What to do?
# identical(graph_orig$from_ctry, graph_nodes$iso3c[graph_orig$from])
# identical(graph_orig$to_ctry, graph_nodes$iso3c[graph_orig$to])
# identical(graph_add$from_ctry, graph_nodes$iso3c[graph_add$from])
# identical(graph_add$to_ctry, graph_nodes$iso3c[graph_add$to])

# Saving
for (name in .c(graph_orig, graph_add, graph_nodes)) {
  sprintf("data/transport_network/%s_MA_100kmh.csv", name) |> 
    fwrite(x = get(name))
}



# STILL NEED TO DO THIS!!
# TODO: Also Adding the information to the RData file
load("data/transport_network/trans_ECA_network.RData")

# Load previous saved graphs
graphs <- sapply(.c(graph_orig, graph_add, graph_nodes), function(name)
  sprintf("data/transport_network/%s_MA_100kmh.csv", name) |> fread())
graphs$graph_orig$add <- FALSE
graphs$graph_add$add <- TRUE

# Joining
nodes %<>% transform(qDF(round(st_coordinates(.), 4))) %>% 
  join(tfm(graphs$graph_nodes, X = round(lon, 4), Y = round(lat, 4)), 
       on = c("X", "Y", "population", "key_city"), drop = "x", overid = 0, how = "inner") %>% select(-X, -Y)
edges %<>% join(graphs$graph_orig, on = c("from", "to"), drop = "x", overid = 2)
add_links %<>% join(graphs$graph_add, on = c("from", "to"), drop = "x", overid = 2)

# Check that network aligns with nodes
allv(st_distance(st_geometry(net, "nodes"), st_geometry(nodes), by_element = TRUE), 0)
identical(st_geometry(net, "edges"), st_geometry(edges))

# Add to network
net %<>% activate("nodes") %>% dplyr::mutate(nodes |> atomic_elem() |> qDF())
net %<>% activate("edges") %>% dplyr::mutate(select(edges, -(from:gravity_dur)) |> atomic_elem() |> qDF())

# Save
TAN_env <- new.env()
load("data/transport_network/trans_ECA_network.RData", envir = TAN_env)
TAN_env$nodes_param <- nodes
TAN_env$edges_param <- edges
TAN_env$add_links_param <- add_links
TAN_env$net_param <- net
save(list = ls(TAN_env), file = "data/transport_network/trans_ECA_network_param_MA_100kmh.RData", envir = TAN_env)



