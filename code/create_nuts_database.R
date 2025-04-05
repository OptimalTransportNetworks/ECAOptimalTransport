library(fastverse)
library(sf)

# Load geometries
NUTS <- geojsonsf::geojson_sf("data/NUTS_RG_01M_2016_3035.geojson") |> 
  janitor::clean_names()
# input = st_crs("EPSG:3035"), wkt = "EPSG:3035"
st_crs(NUTS) <- 3035
NUTS %<>% st_make_valid()

# NUTS |> subset(levl_code == 1) |> mapview::mapview()
# NUTS |> subset(levl_code == 2) |> mapview::mapview()
# NUTS |> subset(levl_code == 3) |> mapview::mapview()

# Load population data
NUTS3_pop <- fread("data/estat_demo_r_pjanaggr3.tsv.gz", header = TRUE, na.strings = ":") %>%
  ftransform(tstrsplit(.[[1]], ",") |> set_names(strsplit(names(.)[1], ",")[[1]])) |> 
  fsubset(sex == "T" & age %in% c("TOTAL", "Y15-64"), -1) |>
  ftransformv(`1990`:`2024`, function(x) if(is.character(x)) as.numeric(tstrsplit(x, " ")[[1]]) else x) |>
  rm_stub("\\TIME_PERIOD", pre = FALSE) |>
  get_vars(varying) |>
  pivot("geo", names = "age", how = "w", transpose = "names", check.dups = TRUE) |>
  colorderv(is_categorical) |>
  janitor::clean_names()

fnobs(NUTS3_pop)

# Joining
NUTS %<>% join(NUTS3_pop, on = c("nuts_id" = "geo")) %>% 
  get_vars(function(x) is.list(x) || varying(x))

fnobs(NUTS)

# NUTS |> subset(levl_code == 3, y15_64_2019) |> mapview::mapview()

# Saving
saveRDS(NUTS, "data/NUTS3_pop.rds")

# Using tmap package, draw a map of the population in 2019
library(tmap)
NUTS_simplified <- rmapshaper::ms_simplify(NUTS) |> st_transform(4326)
europe_bbox <- st_bbox(c(xmin = -12, xmax = 31, ymax = 71, ymin = 35), crs = st_crs(4326))

pdf("figures/population_15_64_NUTS3_2019.pdf", width = 16.5, height = 8)
tm_basemap("Esri.WorldGrayCanvas", zoom = 5) + # CartoDB.Positron
  tm_shape(NUTS_simplified |> subset(levl_code >= 1) |>
             transform(y15_64_2019 = pmax(1, y15_64_2019), 
                       nuts_level = paste("NUTS Level", levl_code)), 
           bbox = europe_bbox) +
  tm_polygons(col = tm_const(), lwd = 0.1, col.legend = tm_legend_hide(),
              fill = "y15_64_2019", 
              fill.scale = tm_scale_continuous_sqrt(values = "brewer.yl_or_rd",
                                               ticks = c(1, 2.5e5, 5e5, 1e6, 2.5e6, 10e6)),
              fill.legend = tm_legend("2019 Population Ages 15-64", position = c("left", "top"), 
                                     stack = "h", frame = FALSE, height = 10, item.width = 0.5, 
                                     text.size = 1.2, title.size = 1.3, title.padding = c(-0.5, 0, 0, 0))) +
  tm_facets_wrap(by = "nuts_level", nrow = 1) +
  tm_layout(frame = FALSE) 
dev.off()

