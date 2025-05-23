library(fastverse)
set_collapse(mask = c("manip", "helper", "special"), nthreads = 4)
fastverse_extend(qs, sf, s2, units, install = TRUE)
source("code/fs_build/helpers.R")
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
# F&S set delta0 so that K = 1000 in each country
get_delta0_ctry <- function(ctry) {
  base <- 0.12 * log(shape$rugg+1) + 0.085 * log(shape$pop_wpop_km2+1)
  ind <- which(shape$iso[routes$from] == ctry & shape$iso[routes$to] == ctry)
  objective <- function(delta0) {
    cost_kmh <- exp(log(delta0) + base)
    cost_kmh <- (cost_kmh[routes$from] + cost_kmh[routes$to]) / 2
    abs(sum(cost_kmh[ind] * routes$time_efficiency[ind]) - 1000)
  }
  c(optimise(objective, c(0, 1e8), tol = 1e-12), list(n = length(ind)))
}

costs <- sapply(unique(shape$iso), get_delta0_ctry, simplify = FALSE) |> 
  rowbind(idcol = "ctry", return = "data.frame")

# Old: 
# shape$cost_km <- exp(log(120e3) + 0.12 * log(shape$rugg+1) + 0.085 * log(shape$pop_wpop_km2+1))
# routes$cost_km <- (shape$cost_km[routes$from] + shape$cost_km[routes$to]) / 2
# routes$cost <- routes$cost_km * routes$sp_distance / 1e6
shape$cost_km <- exp(log(costs$minimum[ckmatch(shape$iso, costs$ctry)]) + 0.12 * log(shape$rugg+1) + 0.085 * log(shape$pop_wpop_km2+1))
routes$cost_km <- (shape$cost_km[routes$from] + shape$cost_km[routes$to]) / 2
routes$cost <- routes$cost_km * routes$time_efficiency
descr(routes$cost_km)
descr(routes$cost)

# Now Getting 120km buffer around each country
countries <- wbstats::wb_cachelist$countries %$% iso3c[region_iso3c %==% "ECS"] %>% 
  c("XKO") %>% setdiff(c("GRL", "ISL", "GBR", "IRL", "NOR", "SWE", "FIN", "FRO", "IMN", "CYP")) %>% 
  intersect(shape$iso)
table(shape$iso)

# No buffers version
for (c in countries) {
  cat(c, " ")
  ind <- unique(shape$iso %==% c)
  if(length(ind) <= 2) next
  cs <- subset(qDT(shape), ind, cell_id, subcell_id, iso, pwx, pwy, 
               predicted_GCP_const_2017_USD, predicted_GCP_const_2017_PPP, pop_cell, national_population,
               cell_GDPC_const_2017_USD, cell_GDPC_const_2017_PPP, is_cell_censored, 
               open_water, rugg, pop_wpop, pop_wpop_km2, cost_km) 

  # Creating product specification
  cs$own_product = cs %in% subset(pop_cell > 200e3,
     .x = largest_within_radius(cs, c("pwx", "pwy"), size = "pop_cell", radius_km = 80))
  print(table(cs$own_product))
  
  if(sum(cs$own_product) > 15) cs[own_product == TRUE, own_product := replace(own_product, -topn(pop_cell, 15), FALSE)]
  cs[own_product == FALSE, product := unattrib(cut(pop_cell, quantile(pop_cell, seq(0, 1, 0.2)), include.lowest = TRUE))]
  cs[own_product == TRUE, product := seq_along(pop_cell) + 5L]
  if(anyNA(cs$product)) stop("Missing product")
  
  fwrite(cs, sprintf("data/grid_network/country_nobuff/%s_nodes.csv", c))
  
  subset(qDT(routes), from %in% ind & to %in% ind, 
         from:to_lat, fx:ty, duration, distance, sp_distance:cost) |> 
    mutate(from = ckmatch(from, ind), to = ckmatch(to, ind)) |> 
    fwrite(sprintf("data/grid_network/country_nobuff/%s_edges.csv", c))
}

# Takes long for Russia
buffers <- sapply(countries, function(x) {
  ind <- unique(shape$iso %==% x)
  cs <- ss(shape, ind)
  dmat <- st_distance(shape, cs) %/=% 1e3
  lapply(mctl(dmat), function(y) which(y < 100)) |> 
    do.call(what = "c") |> unique(sort = TRUE) |>
    setdiff(ind)
}, simplify = FALSE)

# Generate Country Files
for (c in countries) {
  cat(c, " ")
  ind <- shape$iso %==% c
  if(length(ind) <= 2) next
  ind <- unique(c(ind, buffers[[c]]))
  cs <- subset(qDT(shape), ind, cell_id, subcell_id, iso, pwx, pwy, 
         predicted_GCP_const_2017_USD, predicted_GCP_const_2017_PPP, pop_cell, national_population,
         cell_GDPC_const_2017_USD, cell_GDPC_const_2017_PPP, is_cell_censored, 
         open_water, rugg, pop_wpop, pop_wpop_km2, cost_km) |> 
        mutate(is_buff = ind %in% buffers[[c]])
  
  # Creating product specification
  cs$own_product = cs %in% subset(pop_cell > 200e3,
    .x = largest_within_radius(cs, c("pwx", "pwy"), size = "pop_cell", radius_km = 80))
  print(table(cs$own_product))
  
  if(sum(cs$own_product) > 15) {
    .opr <- copy(cs$own_product)
    cs[own_product == TRUE, own_product := replace(own_product, -topn(pop_cell, 15), FALSE)]
    if(cs[, sum(own_product & !is_buff)] < 10) {
      cs[.opr & !is_buff, own_product := replace(own_product, topn(pop_cell, 10), TRUE)]
    }
    rm(.opr)
  }
  cs[own_product == FALSE, product := unattrib(cut(pop_cell, quantile(pop_cell, seq(0, 1, 0.2)), include.lowest = TRUE))]
  cs[own_product == TRUE, product := seq_along(pop_cell) + 5L]
  if(anyNA(cs$product)) stop("Missing product")

  fwrite(cs, sprintf("data/grid_network/country/%s_nodes.csv", c))
  
  subset(qDT(routes), from %in% ind & to %in% ind, 
         from:to_lat, fx:ty,  duration, distance, sp_distance:cost) |> 
  mutate(is_buff = from %in% buffers[[c]] | to %in% buffers[[c]], 
         from = ckmatch(from, ind), to = ckmatch(to, ind)) |> 
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


