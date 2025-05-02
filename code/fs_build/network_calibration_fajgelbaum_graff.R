library(fastverse)
set_collapse(mask = c("manip", "helper", "special"), nthreads = 4)
fastverse_extend(qs, sf, units, install = TRUE)
source("code/helpers.R")
fastverse_conflicts()

shape <- qread("data/grid_network/cells_shape.qs")
routes <- qread("data/grid_network/routes_raw_pw.qs")
segments <- qread("data/grid_network/segments.qs") |> 
  join(routes, on = c("from", "to"), how = "semi")
if(!identical(.subset(routes, .c(from, to)), .subset(segments, .c(from ,to)))) stop("Mismatch")

segments_buff <- st_buffer(segments, as_units(3000, "m"))

# Open Water: Remove certain routes
open_water <- terra::rast("/Users/sebastiankrantz/Documents/Data/Landcover/Consensus_reduced_class_12_open_water.tif")
shape$open_water <- exactextractr::exact_extract(open_water, shape, fun = "mean")
segments$open_water <- routes$open_water <- exactextractr::exact_extract(open_water, segments_buff, fun = "mean")
rm(open_water)
descr(shape$open_water)
descr(routes$open_water)
# mapview::mapview(subset(segments, open_water > 50 & from != 5316)) # 5316 = crossing to Krim
routes %<>% subset(open_water <= 50 | from == 5316)
# mapview::mapview(routes)
routes %<>% ss(-c(23658, 28219, 28220, 27470, 28232, 28717, 25995))
routes %<>% ss(-c(23657)) # In Italy: tiny Island
routes %<>% ss(-c(23656))

# Adding Ruggedness: https://diegopuga.org/data/rugged/
rugg <- terra::rast("/Users/sebastiankrantz/Documents/Data/Ruggedness/tri.txt")
rugg_area <- terra::rast("/Users/sebastiankrantz/Documents/Data/Ruggedness/cellarea.txt")
shape$rugg <- exactextractr::exact_extract(rugg, shape, weights = rugg_area, fun = "weighted_mean")
rm(rugg, rugg_area)

# Adding Population (WorldPop 2020 1km2 global)
pop_wpop <- terra::rast("/Users/sebastiankrantz/Documents/Data/WorldPop/ppp_2020_1km_Aggregated.tif")
shape$pop_wpop <- exactextractr::exact_extract(pop_wpop, shape, fun = "sum")
shape$pop_wpop_km2 <- unattrib(shape$pop_wpop / (st_area(shape) / 1e6))
rm(pop_wpop)

# Computing infrastructure building cost following Collier et. al. (2016)
shape$cost_km <- exp(log(120e3) + 0.12 * log(shape$pop_wpop_km2+1) + 0.085 * log(shape$rugg+1))
routes$cost_km <- (shape$cost_km[routes$from] + shape$cost_km[routes$to]) / 2
routes$cost <- routes$cost_km * routes$sp_distance / 1e6

# Now Getting 120km buffer around each country
countries <- wbstats::wb_cachelist$countries %$% iso3c[region_iso3c %==% "ECS"] %>% 
  c("XKO") %>% setdiff(c("GRL", "ISL", "GBR", "IRL", "NOR", "SWE", "FIN", "FRO", "IMN", "CYP")) %>% 
  intersect(shape$iso)
table(shape$iso)

# Takes long for Russia
buffers <- sapply(countries, function(x) {
  ind <- shape$iso %==% x
  cs <- ss(shape, ind)
  dmat <- st_distance(shape, cs) %/=% 1e3
  lapply(mctl(dmat), function(y) which(y < 100)) |> 
    do.call(what = "c") |> unique(sort = TRUE) |>
    setdiff(ind)
}, simplify = FALSE)

# Generate Country Files
for (c in countries) {
  cat(c, " ")
  ind <- unique(c(shape$iso %==% c, buffers[[c]]))
  subset(qDT(shape), ind, cell_id, subcell_id, iso, pwx, pwy, 
         predicted_GCP_const_2017_USD, predicted_GCP_const_2017_PPP, pop_cell, national_population,
         cell_GDPC_const_2017_USD, cell_GDPC_const_2017_PPP, is_cell_censored, 
         open_water, rugg, pop_wpop, pop_wpop_km2, cost_km) |> 
    fwrite(sprintf("data/grid_network/country/%s_nodes.csv", c))
  subset(qDT(routes), from %in% ind | to %in% ind, 
               from:to_lat, fx:ty,  duration, distance, sp_distance:cost) |> 
    fwrite(sprintf("data/grid_network/country/%s_edges.csv", c))
}

# Save Overall
select(qDT(shape), cell_id, subcell_id, iso, pwx, pwy, 
       predicted_GCP_const_2017_USD, predicted_GCP_const_2017_PPP, pop_cell, national_population,
       cell_GDPC_const_2017_USD, cell_GDPC_const_2017_PPP, is_cell_censored, 
       open_water, rugg, pop_wpop, pop_wpop_km2, cost_km) |> 
  fwrite("data/grid_network/ECA_nodes.csv")

select(qDT(routes), from:to_lat, fx:ty,  duration, distance, sp_distance:cost) |> 
  fwrite("data/grid_network/ECA_edges.csv")


