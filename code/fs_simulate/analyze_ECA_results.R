#####################################################################
# Transport Network: Analyze Optimal Trans-ECA GE Investments
#####################################################################

library(fastverse)
set_collapse(mask = c("manip", "helper", "special"))
fastverse_extend(qs, sf, units, sfnetworks, tmap, install = TRUE)
# source("code/helpers.R")
fastverse_conflicts()

res_name <- "ug_29g_25b_fixed_irs_cgc_sigma3.8_rho2" # "20g_2327m_fixed_cgc_sigma3.8_rho0_julia_google" # '22g_add_10b_fixed_duality_sigma3'

results <- list(
  nodes = fread(sprintf("results/transport_network/GE/nodes_results_%s.csv", res_name)) |> tfmv(c(lon, lat), round, 5),
  edges = fread(sprintf("results/transport_network/GE/edges_results_%s.csv", res_name))
)

network <- new.env()
load("data/transport_network/trans_ECA_network.RData", envir = network) 
nodes <- network$nodes |> tfm(qDT(st_coordinates(geometry)) |> set_names(c("lon", "lat"))) |> tfmv(c(lon, lat), round, 5)
edges <- qread("data/transport_network/edges_real_simplified.qs")
tfm(edges) <- qDT(network$edges) |> select(-geometry)

largest <- results$nodes %>% subset(unclass(product) > 4L) %$% set_names(product, city_ascii) %>% sort() %>% names()
attr(results$nodes$product, "levels") <- c("Small City", "City > 100K", "City > 250K", "City > 1M", largest)
class(results$nodes$product) <- "factor"

if(!res_name %ilike% "ug_") {
  edges <- rbind(
    select(network$add_links, from, to) |> tfm(add = TRUE),
    select(edges, from, to) |> tfm(add = FALSE))
  results$edges %<>% join(x = edges, on = c("from", "to", "add"), how = "inner", drop = "x", overid = 2)
  # results$edges$add[between(results$edges$Ijk_orig, 14.2948, 14.2949)] <- TRUE
} else {
  results$edges %<>% join(x = edges, on = c("from", "to"), how = "inner", drop = "x")
}
results$nodes %<>% join(x = nodes, on = c("lon", "lat"), how = "inner", drop = "x")
results$edges %<>% mutate(distance = distance, # / 1000,
                          cost_per_km = total_cost / distance)

# Statistics on the upgrade extent -----------------------------------------------------------------

# Cost of all possible work (should be 13.5 billion in millions)
results$edges |> with(sum(distance*cost_per_km))
# Budget spent (should give 1 or 2 billion in millions)
results$edges |> with(sum(replace_inf(pmax((Ijk-Ijk_orig)/(pmax(Ijk_orig, 100)-Ijk_orig), 0), 0)*distance*cost_per_km))
# Amount spent on different types of work (1 = new construction, 0 = upgrade)
# results$edges |> with(fsum(replace_inf(pmax((Ijk-Ijk_orig)/(pmax(Ijk_orig, 100)-Ijk_orig), 0), 0)*distance*cost_per_km, add)) # |> proportions()

# Road km built/upgraded
results$edges |> with(sum(replace_inf(pmax((Ijk-Ijk_orig)/(pmax(Ijk_orig, 100)-Ijk_orig), 0), 0)*distance))
# Km on different types of work
results$edges |> with(fsum(replace_inf(pmax((Ijk-Ijk_orig)/(pmax(Ijk_orig, 100)-Ijk_orig), 0), 0)*distance, add)) # |> proportions()

# Number of roads worked on (TRUE) (extensive margin)
results$edges |> with(Ijk-Ijk_orig > 1) |> table()  # |> proportions()
results$edges |> with(table(Ijk-Ijk_orig > 1, add)) # proportions() |> addmargins() # By type

# Work Intensity (intensive margin = km/h added)
results$edges |> with(descr((Ijk-Ijk_orig)[Ijk-Ijk_orig > 1]))
results$edges |> subset(Ijk-Ijk_orig > 1) |> with(descr(Ijk-Ijk_orig, add)) # By type


# Statistics on the economic gains -----------------------------------------------------------------

if(res_name %ilike% "rho2") results$nodes %<>% mutate(uj = (uj*(-1))^(-1), uj_orig = (uj_orig*(-1))^(-1))

# Global Welfare Gains (Ratio)
results$nodes |> with(sum(uj * Lj) / sum(uj_orig * Lj_orig))
# Local Welfare Gains (Ratio)
results$nodes |> with(descr(uj / uj_orig))
# Correlates of Local Welfare Gains (Ratio)
qDF(results$nodes) |> mutate(ugain = uj / uj_orig) |>
  select(ugain, population, GDPC, GDP) |> pwcor()

# Consumption Gains
results$nodes |> with(sum(Cj) / sum(Cj_orig)) # Global
results$nodes |> with(descr(Cj / Cj_orig))    # Local

# Consumption Percent by City Type: May need to change Dj -> Cj if cgc is off (all_routes/dual solutions)
# Overall
qDT(results$nodes) |> group_by(product) |> gvr("^Dj_") |> fsum() |> 
  tfmv(-1, fsum, TRA = "%", apply = FALSE) |> qM(1) %>% {set_names(diag(.), rownames(.))}
# At the city level + median aggregation 
qDT(results$nodes) |> group_by(product) |> gvr("^Dj_") %>%
  dapply(`/`, psum(.)) |> fmedian() |> qM(1) %>% {set_names(diag(.), rownames(.))}
# Per capita
qDT(results$nodes) |> group_by(product) |> gvr("^Dj_|pop") |> 
  tfmv(-population, `/`, population+1) |>
  fmedian() |> tfmv(-(1:2), fsum, TRA = "%", apply = FALSE) |> select(-population) |>
  qM(1) %>% {set_names(diag(.), rownames(.))}


# Plots of Final Network and Optimal Investments ---------------------------------------------------

results$nodes %<>% mutate(prod2 = droplevels(fifelse(unclass(product) > 4L, NA, product)))

# tmap_options(raster.max_cells = 1e6)

# Final Network
pdf(sprintf("figures/GE/trans_ECA_network_GE_%s.pdf", res_name), width = 10, height = 4.2)
tm_basemap("CartoDB.Positron", zoom = 5) +
  tm_shape(results$edges) +
  tm_lines(col = "Ijk", 
           col.scale = tm_scale_continuous(7, values = "inferno", limits = c(0, 110)),
           col.legend = tm_legend("km/h", position = c("left", "top"), frame = FALSE, 
                                  text.size = 1, title.size = 1.2, margins = c(0, -0.5, 0, 0),
                                  title.padding = c(0, 0, -0.5, 0), 
                                  item.space = 0, item.height = 1.5, item.width = 0.5), lwd = 1.5) +
  tm_shape(results$nodes) + 
  tm_dots(fill = "prod2", size = 0.1, 
          fill.scale = tm_scale_categorical(values = "turbo", value.na = "purple3", label.na = "Large City"),
          fill.legend = tm_legend("Product", position = c("right", "bottom"), 
                                  frame = FALSE, text.size = 0.75, title.size = 1, 
                                  margins = c(-0.5, 0, 0, -0.5))) +
  tm_shape(subset(results$nodes, is.na(prod2))) + tm_dots(size = 0.25, fill = "purple3") +
  tm_layout(frame = FALSE)
dev.off()

# Difference (Intensive Margin)
pdf(sprintf("figures/GE/trans_ECA_network_GE_%s_diff.pdf", res_name), width = 10, height = 4.2)
tm_basemap("CartoDB.Positron", zoom = 5) +
  tm_shape(mutate(results$edges, diff = pmin(pmax(Ijk - Ijk_orig, 0), 100))) +
  tm_lines(col = "diff", 
           col.scale = tm_scale_continuous(values = "brewer.yl_or_rd"),
           col.legend = tm_legend(expression(Delta~"km/h"), position = c("left", "top"), frame = FALSE, 
                                  text.size = 1, title.size = 1.2, margins = c(0, -0.5, 0, 0),
                                  title.padding = c(0, 0, -0.5, 0), 
                                  item.space = 0, item.height = 1.5, item.width = 0.5), lwd = 1.5) +
  tm_shape(results$nodes) + 
  tm_dots(fill = "prod2", size = 0.1, 
          fill.scale = tm_scale_categorical(values = "turbo", value.na = "purple3", label.na = "Large City"),
          fill.legend = tm_legend("Product", position = c("right", "bottom"), 
                                  frame = FALSE, text.size = 0.75, title.size = 1, 
                                  margins = c(-0.5, 0, 0, -0.5))) +
  tm_shape(subset(results$nodes, is.na(prod2))) + tm_dots(size = 0.25, fill = "purple3") +
  tm_layout(frame = FALSE)
dev.off()

# Upgrade Percent (Extensive Margin)
pdf(sprintf("figures/GE/trans_ECA_network_GE_%s_perc_ug.pdf", res_name), width = 10, height = 4.2)
tm_basemap("CartoDB.Positron", zoom = 5) +
  tm_shape(mutate(results$edges, perc_ug = pmin(pmax((Ijk - Ijk_orig)/(100 - Ijk_orig)*100, 0), 100))) +
  tm_lines(col = "perc_ug", 
           col.scale = tm_scale_continuous(ticks = seq(0, 100, 20), values = "brewer.yl_or_rd"), 
           col.legend = tm_legend("% Completed", position = c("left", "top"), frame = FALSE, 
                                  text.size = 1, title.size = 1.2, margins = c(0, -0.5, 0, 0),
                                  title.padding = c(0, 0, -0.5, 0), 
                                  item.space = 0, item.height = 1.5, item.width = 0.5), lwd = 1.5) +
  tm_shape(results$nodes) + 
  tm_dots(fill = "prod2", size = 0.1, 
          fill.scale = tm_scale_categorical(values = "turbo", value.na = "purple3", label.na = "Large City"),
          fill.legend = tm_legend("Product", position = c("right", "bottom"), 
                                  frame = FALSE, text.size = 0.75, title.size = 1, 
                                  margins = c(-0.5, 0, 0, -0.5))) +
  tm_shape(subset(results$nodes, is.na(prod2))) + tm_dots(size = 0.25, fill = "purple3") +
  tm_layout(frame = FALSE)
dev.off()

# Flow of Goods: 4 Cities
pdf(sprintf("figures/GE/trans_ECA_network_GE_%s_good_flows_4_city.pdf", res_name), width = 11, height = 5)
tm_basemap("CartoDB.Positron", zoom = 5) +
  tm_shape(results$edges |> gvr("Qjk_") |> pivot("geometry") |> 
             mutate(variable = set_attr(variable, "levels", tools::toTitleCase(levels(results$nodes$product)))) |>
             subset(variable %ilike% "Moscow|Istanbul|Paris|Tashkent")) +
  tm_facets_wrap("variable", ncol = 2) + 
  tm_lines(col = "value", 
           col.scale = tm_scale_continuous(8, trans = "log1p", values = "yl_or_rd"),
           col.legend = tm_legend("Flow", position = tm_pos_in(0, 1), height = 10.5, 
                                  title.size = 0.75, text.size = 0.55, frame = FALSE), lwd = 1.5) +
  tm_layout(frame = FALSE, 
            panel.label.bg.color = "white",
            panel.label.frame = FALSE) 
dev.off()

# Local Welfare (Utility Per Worker) Gains
pdf(sprintf("figures/GE/trans_ECA_network_GE_%s_upw_gain.pdf", res_name), width = 10, height = 4.2)
tm_basemap("CartoDB.Positron", zoom = 5) +
  tm_shape(results$edges) +
  tm_lines(col = "grey", lwd = 1.5) +
  tm_shape(mutate(results$nodes, ugain = (uj / uj_orig - 1) * 100) |>
             st_as_sf(coords = .c(lon, lat), crs = 4326)) +
  tm_dots(fill = "ugain", 
          fill.scale = tm_scale_intervals(7, values = "turbo"),
          fill.legend = tm_legend("Welfare Gain (%)", position = c("left", "top"), frame = FALSE, 
                                  text.size = 1, title.size = 1.2, margins = c(0, -0.5, 0, 0),
                                  title.padding = c(0, 0, -0.5, 0), 
                                  item.space = 0, item.height = 1.5, item.width = 0.5), size = 0.2) +
  tm_layout(frame = FALSE) 
dev.off()






# Flows of Goods -----------------------------------------------------------------------------------

tmap_options(raster.max.cells = 1e7)

# Flow of Goods: All Cities
pdf(sprintf("figures/GE/trans_ECA_network_GE_%s_good_flows.pdf", res_name), width = 10, height = 11)
tm_basemap("CartoDB.Positron", zoom = 5) +
tm_shape(results$edges |> gvr("Qjk_") |> pivot("geometry") |> 
           mutate(variable = set_attr(variable, "levels", tools::toTitleCase(levels(results$nodes$product))))) +
  tm_facets_wrap("variable", ncols = 5) + 
  tm_lines(col = "value", 
           col.scale = tm_scale_continuous(8, trans = "log1p", values = "yl_or_rd"),
           col.legend = tm_legend("Flow", position = tm_pos_in(0.73, 0.47), height = 10, 
                                  title.size = 0.75, text.size = 0.55, frame = FALSE), lwd = 1.5) +
  tm_layout(frame = FALSE, 
            panel.label.bg.color = "white",
            panel.label.frame = FALSE) 
dev.off()

# Flow of Goods: 4 Cities
pdf(sprintf("figures/GE/trans_ECA_network_GE_%s_good_flows_4_city.pdf", res_name), width = 11, height = 5)
tm_basemap("CartoDB.Positron", zoom = 5) +
  tm_shape(results$edges |> gvr("Qjk_") |> pivot("geometry") |> 
             mutate(variable = set_attr(variable, "levels", tools::toTitleCase(levels(results$nodes$product)))) |>
             subset(variable %ilike% "Moscow|Istanbul|Paris|Tashkent")) +
  tm_facets_wrap("variable", ncol = 2) + 
  tm_lines(col = "value", 
           col.scale = tm_scale_continuous(8, trans = "log1p", values = "yl_or_rd"),
           col.legend = tm_legend("Flow", position = tm_pos_in(0, 1), height = 10.5, 
                                  title.size = 0.75, text.size = 0.55, frame = FALSE), lwd = 1.5) +
  tm_layout(frame = FALSE, 
            panel.label.bg.color = "white",
            panel.label.frame = FALSE) 
dev.off()


# Impact of Frictions on Optimal Investments -------------------------------------------------------

tmap_options(raster.max.cells = 1e6)

res_name <- "add_20g_3400m_fixed_cgc_sigma3.8_rho0_julia_MACR_100kmh_google"
edges_res <- list(NoFR = fread(sprintf("results/transport_network/GE/MACR/edges_results_%s.csv", res_name)),
                  FR = fread(sprintf("results/transport_network/GE/MACR/edges_results_%s.csv", sub("_julia", "_bc_julia", res_name))))
edges_res$FR$distance = edges_res$NoFR$distance

# Road km built/upgraded
edges_res |> lapply(with, sum(replace_inf(pmax((Ijk-Ijk_orig)/(pmax(Ijk_orig, 100)-Ijk_orig), 0), 0)*distance)) |>
  mutate(diff = FR - NoFR, perc_diff = diff / NoFR * 100) |> qDF()
# Km on different types of work
if(res_name %ilike% "add") {
  edges_res$NoFR$add[between(edges_res$NoFR$Ijk_orig, 14.2948, 14.2949)] <- TRUE
  edges_res$FR$add[between(edges_res$FR$Ijk_orig, 14.2948, 14.2949)] <- TRUE
}
edges_res |> lapply(with, fsum(replace_inf(pmax((Ijk-Ijk_orig)/(pmax(Ijk_orig, 100)-Ijk_orig), 0), 0)*distance, add)) |> 
  mutate(diff = FR - NoFR, perc_diff = diff / NoFR * 100) |> qM()

# Changes in Welfare
nodes_res <- list(NoFR = fread(sprintf("results/transport_network/GE/MACR/nodes_results_%s.csv", res_name)),
                  FR = fread(sprintf("results/transport_network/GE/MACR/nodes_results_%s.csv", sub("_julia", "_bc_julia", res_name))))
if(res_name %ilike% "rho2") nodes_res %<>% lapply(mutate, uj = (uj*(-1))^(-1), uj_orig = (uj_orig*(-1))^(-1))
nodes_res |> lapply(with, (sum(uj * Lj) / sum(uj_orig * Lj_orig) - 1) * 100) |>
  mutate(diff = FR - NoFR, perc_diff = diff / NoFR * 100) |> qDF()

# Computing Percent Upgraded Difference
edges_res %<>% lapply(select, from, to, Ijk, Ijk_orig) %>% 
  rowbind(idcol = "data") %>% 
  mutate(perc_ug = pmin(pmax((Ijk - Ijk_orig)/(100 - Ijk_orig)*100, 0), 100)) %>% 
  pivot(c("from", "to"), "perc_ug", "data", how = "w") %>% 
  mutate(perc_ug_diff = replace_outliers(FR - NoFR, c(-100, 100), "clip"))
edges_res %<>% join(x = edges, on = c("from", "to"))

descr(edges_res$perc_ug_diff)
replace_na(edges_res$perc_ug_diff, set = TRUE)

# Different in % Upgraded (Extensive Margin)
pdf(sprintf("figures/GE/trans_ECA_network_GE_%s_Ijk_bc_perc_ug_diff.pdf", res_name), width = 10, height = 4.2)
tm_basemap("Esri.WorldGrayCanvas", zoom = 5) + 
  tm_shape(edges_res) +
  tm_lines(col = "perc_ug_diff", 
           col.scale = tm_scale_continuous(ticks = c(-100, -50, -25, 0, 25, 50, 100), 
                                           limits = c(-100, 100), midpoint = 0, values = "-classic_red_blue"), # "-spectral"
           col.legend = tm_legend("Diff in % UG", position = c("right", "bottom"), frame = FALSE,
                                  text.size = 1.25, title.size = 1.7), lwd = 2) +
  tm_layout(frame = FALSE, 
            panel.label.bg.color = "white",
            panel.label.frame = FALSE) 
dev.off()


# Impact of Frictions on Trade Flows ---------------------------------------------------------------

tmap_options(raster.max.cells = 1e7)

res_name <- "add_20g_1b_fixed_cgc_sigma3.8_rho0_julia_google"
edges_res <- list(NoFR = fread(sprintf("results/transport_network/GE/edges_results_%s.csv", res_name)),
                  FR = fread(sprintf("results/transport_network/GE/edges_results_%s.csv", sub("_julia", "_bc_julia", res_name))))
edges_res %<>% lapply(gvr, "^from$|^to$|^Qjk_") # %>% rowbind(idcol = "data")
edges_res$Ratio <- edges_res$FR %>% tfm(slt(., Qjk_1:Qjk_22) %c/% slt(join(slt(edges_res$FR, from, to), edges_res$NoFR), Qjk_1:Qjk_22) %>% 
                                          replace_outliers(c(0, 100), "clip"))
edges_res %<>% lapply(join, x = edges, on = c("from", "to"))

descr(atomic_elem(edges_res$Ratio))
edges_res %>% lapply(. %>% atomic_elem() %>% num_vars() %>% fsum())

pdf(sprintf("figures/GE/trans_ECA_network_GE_%s_good_flows_bc_ratio.pdf", res_name), width = 10, height = 10)
tm_basemap("CartoDB.Positron", zoom = 5) +
  tm_shape(edges_res$Ratio |> gvr("Qjk_") |> pivot("geometry") |> 
           mutate(variable = set_attr(variable, "levels", stringi::stri_trans_general(levels(results$nodes$product), "latin-ascii"))) |> na_omit()) +
  tm_facets_wrap("variable", nrows = 5) + 
  tm_lines(col = "value", 
           col.scale = tm_scale_intervals(breaks = c(0, 0.25, 0.5, 0.75, 1, 1.33, 2, 4, Inf), 
                                          midpoint = 1, values = "-rd_bu"),
           col.legend = tm_legend("Flows Ratio", position = tm_pos_in(0.02, 0.54), height = 10, 
                                  title.size = 0.75, text.size = 0.55, frame = FALSE), lwd = 2) +
  tm_layout(frame = FALSE, 
            panel.label.bg.color = "white",
            panel.label.frame = FALSE) 
dev.off()


# Saving All results ---------------------------------------------------------------


files <- list.files("results/transport_network/GE") |> 
         grep(pattern = "google", value = TRUE) |>
         grep(pattern = "edges", value = TRUE) |>
         grep(pattern = "add", invert = TRUE, value = TRUE)

results <- sapply(files, function(res_name) {
  fread(sprintf("results/transport_network/GE/%s", res_name)) |>
    select(from, from_ctry, to, to_ctry, upgrade_cat, ug_cost_km, Ijk_orig, Ijk)
  }, simplify = FALSE) |> rowbind(idcol = "data") |>
  mutate(data = gsub("_julia_google|_fixed|_cgc|_20g|edges_results|.csv|_|\\.", "", data), 
         data = sub("sigma38", "", data), 
         data = paste0("Ijk_", sub("rho", "r", data))) |>
  pivot(2:8, "Ijk", "data", how = "w", sort = TRUE, check = TRUE)

edges %<>% join(results)

# pmin(pmax((Ijk - Ijk_orig)/(100 - Ijk_orig)*100, 0), 100))

# Add links: 
files <- list.files("results/transport_network/GE") |> 
  grep(pattern = "google", value = TRUE) |>
  grep(pattern = "edges", value = TRUE) |>
  grep(pattern = "add", value = TRUE)

results <- sapply(files, function(res_name) {
  fread(sprintf("results/transport_network/GE/%s", res_name)) |> 
    subset(add, from, to, cost_per_km, Ijk_orig, Ijk)
}, simplify = FALSE) |> rowbind(idcol = "data") |>
  mutate(data = gsub("_julia_google|_fixed|_cgc|_20g|edges_results|.csv|_|\\.", "", data), 
         data = sub("sigma38", "", data), 
         data = paste0("Ijk_", sub("rho", "r", data))) |>
  pivot(2:5, "Ijk", "data", how = "w", sort = TRUE, check = TRUE)

network$add_links %<>% join(results)
