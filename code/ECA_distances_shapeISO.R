library(fastverse)
fastverse_extend(sf, haven)

ECA_centroids <- fread("data/ECA_shp_new/ECA_centroids.csv") |> rm_stub("_") |> 
  fsubset(is.finite(CX) & is.finite(CY)) |>
  collap( ~ countrycode + shapeISO)
# -> Matias fixed Kazhakstan, Armenia, and Georgia

ECA_centroids |> 
  fselect(shapeISO, shapeISO_nm, countrycode, lon = CX, lat = CY) |>
  fwrite("data/ECA_centroids.csv")

ECA_centroids %<>% st_as_sf(coords = c("CX", "CY"), crs = 4326)

# plot(ECA_centroids[, "shapeISO"])
# mapview::mapview(ECA_centroids)

# Now Computing Distance Matrix
library(osrm)
split_large_dist_matrix <- function(data, chunk_size = 100, verbose = FALSE) {
  n = nrow(data)
  res_list = list()
  
  # Loop over each chunk to compute the pairwise distances and travel times
  count = 0
  for (i in seq(1, n, by = chunk_size)) {
    for (j in seq(1, n, by = chunk_size)) {
      # Define the row indices for the current chunks
      rows_i = i:min(i + chunk_size - 1, n)
      rows_j = j:min(j + chunk_size - 1, n)
      
      # Extract the data for the current chunks
      ds_i = data[rows_i, ]
      ds_j = data[rows_j, ]
      
      if(verbose) {
        count = count + 1L
        cat(count," ")
      }
      
      # Perform the API call for the current chunks
      r_ij = osrmTable(src = ds_i, dst = ds_j, measure = c('duration', 'distance'))
      
      # Store the result in a list for later combination
      res_list[[paste(i, j, sep = "_")]] = r_ij
    }
  }
  
  # Combine the results from the list into one large matrix for durations and distances
  res_sources = matrix(NA, n, 2)
  res_destinations = matrix(NA, n, 2)
  res_durations = matrix(NA, n, n)
  res_distances = matrix(NA, n, n)
  for (i in seq(1, n, by = chunk_size)) {
    for (j in seq(1, n, by = chunk_size)) {
      rows_i = i:min(i + chunk_size - 1, n)
      rows_j = j:min(j + chunk_size - 1, n)
      
      # Retrieve the result from the list
      r_ij = res_list[[paste(i, j, sep = "_")]]
      
      # Place the result into the corresponding location in the matrix
      res_sources[rows_i, ] = qM(r_ij$sources)
      res_destinations[rows_j, ] = qM(r_ij$destinations)
      res_durations[rows_i, rows_j] = r_ij$durations
      res_distances[rows_i, rows_j] = r_ij$distances
    }
  }
  
  # Create a result list to return
  res = list(
    sources = qDF(copyAttrib(mctl(res_sources), data)),
    destinations = qDF(copyAttrib(mctl(res_destinations), data)),
    durations = res_durations,
    distances = res_distances
  )
  
  rn = rownames(res$sources)
  if(length(rn) && suppressWarnings(!identical(as.integer(rn), seq_along(rn)))) {
    dimnames(res$durations) <- dimnames(res$distances) <- list(rn, rn)
  }
  
  return(res)
}

centroids <- ECA_centroids |> 
  ftransform(st_coordinates(geometry) |> mctl() |> set_names(c("lon", "lat"))) |>
  unclass() |> fselect(shapeISO, lon, lat) |> qM(1) |> qDF()
  
dist <- split_large_dist_matrix(centroids)
anyNA(dist$distances)
anyNA(dist$durations)
dist$destinations <- NULL
dist$centroids <- centroids
setrename(dist, sources = starts)

# Road Transport Cost Estimation
# ChatGPT suggests: Cost [€] = 1.05 × (distance in km) + 0.60 × (travel time in minutes)
# Iimi (2023) suggests: log(U.S. cents per ton-km) = 4.650 - 0.395 × log(speed in km/h) - 0.064 × log(distance in km) + 0.024 × crossborder dummy
# Also sector-specific results and different equations for domestic/international shipments...

# Compute travel cost
dist$cents_per_ton_km <- exp(4.650 - 0.395 * log(dist$distances / dist$durations * 60 / 1000) - 
                               0.064 * log(dist$distances / 1000) + 
                               0.024 * outer(ECA_centroids$countrycode, ECA_centroids$countrycode, "!="))
diag(dist$cents_per_ton_km) <- 0
dist$dollars_per_ton <- dist$cents_per_ton_km * dist$distances / 10000

saveRDS(dist, "data/ECA_centroids_distances_and_costs.rds")

diag(dist$cents_per_ton_km) <- NA
descr(vec(dist$cents_per_ton_km))

diag(dist$dollars_per_ton) <- NA
descr(vec(dist$dollars_per_ton))

result <- dist |>
  atomic_elem() |>
  unlist2d("variable", "from") |>
  pivot("from", names = list(from = "variable", to = "to"), how = "r") |>
  fsubset(from != to) |>
  fmutate(distances = distances / 1000) |>
  frename(shapeISO_o = from,
          shapeISO_d = to,
          distance_km = distances,
          duration_min = durations)

fnobs(result)

result |> fwrite("data/ECA_centroids_distances_and_costs.csv")

# Test joining to database

ECA <- read_dta("data/ECA_database_shapeISO.dta")
ECA |> join(result, on = c("shapeISO" = "shapeISO_o")) |> invisible() 
ECA |> join(result, on = c("shapeISO" = "shapeISO_d")) |> invisible() 


