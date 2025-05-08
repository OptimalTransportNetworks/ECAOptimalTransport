library(fastverse)
fastverse_extend(sf, haven)

ECA_centroids <- fread("data/ECA_centroids.csv")
# ctry_dist <- read_dta("/Users/sebastiankrantz/Documents/Data/geodist/dist_cepii.dta")
ctry <- read_dta("/Users/sebastiankrantz/Documents/Data/geodist/geo_cepii.dta")
cities <- fread("/Users/sebastiankrantz/Documents/Data/cities/simplemaps_worldcities_basicv1.77/worldcities.csv")
setv(cities$iso3, "XKS", "XKX")
all(ECA_centroids$countrycode %in% cities$iso3)
cities %<>% fsubset(iso3 %in% ECA_centroids$countrycode)

centroids <- ECA_centroids |> st_as_sf(coords = c("lon", "lat"), crs = 4326)
cities %<>% st_as_sf(coords = c("lng", "lat"), crs = 4326)
dmat <- st_distance(cities, centroids)
dmat[outer(cities$iso3, centroids$countrycode, "!=")] <- Inf
cities$shapeISO <- ECA_centroids$shapeISO[dapply(dmat, which.min, MARGIN = 1)]
rm(dmat, centroids)

cities %<>% qDF() %>% collap(population ~ shapeISO, fsum)
ECA_centroids %<>% join(cities, on = "shapeISO")
descr(ECA_centroids$population)
ECA_centroids[is.na(population) | population < 1000, population := 1000]
descr(ECA_centroids$population)
rm(cities)

ECA_centroids |> fwrite("data/ECA_centroids.csv")

# Now fetching trade data
setdiff(africamonitor::am_countries_wld$ISO3, ctry$iso3)
setdiff(ctry$iso3, africamonitor::am_countries_wld$ISO3)
recode_char(ctry$iso3, ZAR = "COD", ROM = "ROU", SCG = "SRB", YUG = "SRB", set = TRUE) 
ctry %<>% fsubset(iso3 %in% africamonitor::am_countries_wld$ISO3)

IMF2code <- africamonitor::am_countries_wld$ISO2[africamonitor::am_countries_wld$ISO3 %in% ctry$iso3] |> 
  copyv("EH", "ET") |> c("XK", "TW")

DOT <- rdbnomics::rdb("IMF", "DOT", dimensions = list(
  FREQ = "A",
  INDICATOR = c("TMG_CIF_USD", "TXG_FOB_USD"),
  REF_AREA = IMF2code,
  COUNTERPART_AREA = IMF2code
))

DOT <- DOT |> 
  pivot(ids = c("REF_AREA", "Reference Area", "COUNTERPART_AREA", "Counterpart Reference Area", "period"), 
        values = "value", 
        names = "INDICATOR",
        labels = "Indicator",
        how = "wider", 
        check.dups = TRUE) |> 
  janitor::clean_names() |> 
  ftransform(year = year(period), period = NULL) |> 
  fsubset(year >= 2014 & year <= 2024, -year) |> 
  collapv(is_categorical, cols = is.double) |> 
  replace_na(cols = is.double, set = TRUE) |> 
  fcompute( 
    iso_o = countrycode::countrycode(ref_area, "iso2c", "iso3c") |> replace_na("XKX"), 
    country_o = reference_area,
    iso_d = countrycode::countrycode(counterpart_area, "iso2c", "iso3c") |> replace_na("XKX"),
    country_d = counterpart_reference_area,
    ex_fob_usd = txg_fob_usd,
    im_cif_usd = tmg_cif_usd
  )

DOT |> fwrite("data/DOTS.csv")

# Joining to regions
ECA_ctry_centroids <- rowbind(ECA_centroids, 
      fsubset(ctry, iso3 %!in% ECA_centroids$countrycode, shapeISO = iso3, shapeISO_nm = country, lon, lat) |> 
      fmutate(countrycode = shapeISO, population = 1))

dmat <- ECA_ctry_centroids |> st_as_sf(coords = c("lon", "lat"), crs = 4326) |> st_distance()
dimnames(dmat) <- list(ECA_ctry_centroids$shapeISO, ECA_ctry_centroids$shapeISO)
diag(dmat) <- NA

ECA_ctry_trade <- qDT(dmat, "shapeISO_o") |> 
  pivot("shapeISO_o", names = list("shapeISO_d", "distance"), na.rm = TRUE)

ECA_ctry_trade %<>% join(slt(ECA_ctry_centroids, shapeISO, iso3_o = countrycode, pop_o = population), 
                         on = c("shapeISO_o" = "shapeISO")) %>%
                    join(slt(ECA_ctry_centroids, shapeISO, iso3_d = countrycode, pop_d = population), 
                         on = c("shapeISO_d" = "shapeISO")) 

# Joining trade to distance data
setdiff(ctry$iso3, DOT$iso_d)
setdiff(DOT$iso_d, ctry$iso3)
setdiff(ECA_ctry_trade$iso3_d, DOT$iso_d)
setdiff(DOT$iso_d, ECA_ctry_trade$iso3_d)

ECA_ctry_trade %<>% join(slt(DOT, iso_o, iso_d, ex_fob_usd, im_cif_usd), 
                         on = c("iso3_o" = "iso_o", "iso3_d" = "iso_d"), how = "inner")

# Scaling
ECA_ctry_trade %<>% fmutate(
  gravity_share = fsum(pop_o * pop_d / distance, list(iso3_o, iso3_d), TRA = "/"),
  im_cif_usd_gr = im_cif_usd * gravity_share,
  ex_fob_usd_gr = ex_fob_usd * gravity_share
)

# Saving
ECA_ctry_trade |> fwrite("data/ECA_centroids_world_trade_gravity.csv")

