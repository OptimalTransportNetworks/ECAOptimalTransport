library(fastverse)
fastverse_extend(sf, haven, blackmarbler)

# Load DTA
ECA_NUTS <- read_dta("data/ECA_database_country and nuts_apc.dta") |> funique()
# ECA_NUTS |> write_dta("data/ECA_database_country and nuts_apc.dta")
fndistinct(ECA_NUTS)
ECA_NUTS |> fselect(year, nuts2id, nuts3id, nuts2idd, nuts3idd) |> fnunique()
# ECA_NUTS |> gvr("year|nuts") |> fnunique()

NUTS <- readRDS("data/NUTS3_pop.rds")
setdiff(ECA_NUTS$nuts3id, NUTS$nuts_id)
setdiff(NUTS$nuts_id, ECA_NUTS$nuts3id)

# Fetch nightlights for NUTS1 regions
NUTS1 <- NUTS |> subset(levl_code == 1) |> st_transform(4326)
#### Define NASA bearer token
bearer <- "eyJ0eXAiOiJKV1QiLCJvcmlnaW4iOiJFYXJ0aGRhdGEgTG9naW4iLCJzaWciOiJlZGxqd3RwdWJrZXlfb3BzIiwiYWxnIjoiUlMyNTYifQ.eyJ0eXBlIjoiVXNlciIsInVpZCI6InNlYmtyYW50eiIsImV4cCI6MTc0ODMwNzE5NCwiaWF0IjoxNzQzMTIzMTk0LCJpc3MiOiJodHRwczovL3Vycy5lYXJ0aGRhdGEubmFzYS5nb3YiLCJpZGVudGl0eV9wcm92aWRlciI6ImVkbF9vcHMiLCJhY3IiOiJlZGwiLCJhc3N1cmFuY2VfbGV2ZWwiOjN9.CTnC0Do3llf-wt_hba-OkGJBhXzKRshWBFpK4cx3xtJaLauLS7SiAeUrMR1WqAbq9h4R30E8PUCVnnHRMZd5uf5D6VTXxhhqEVT-4yUt33Ji1Q6pzHUwIb30uqKMkggxb_n4gL-8s3oot4b6PwutiU4bgsZX_YALakTeLHnLg8Ry9r9zVFmOJeorYgCpWRsx0XdqT-_YNmYA8MrK9RNZ3FmFh2Hs3_xN_KPBgMGEb6jpNciBQiir1cjvFQuBPWGP-Eix1PfmBN-DJuKXdp35tKGduJtUTpKxJ-n_ZJLCgJZ5Jw63UewHJduu0kziPy3IYBO8ufMkufv3vv_uSyvK7A"
NL21 <- vector("list", length = nrow(NUTS1))
for (i in seq_along(NL21)[-(1:124)]) {
  print(NUTS1$nuts_name[i])
  NL21[[i]] <- bm_raster(roi_sf = NUTS1[i,1],
                         product_id = "VNP46A4",
                         date = "2021",
                         bearer = bearer) |>
               as.data.frame(xy = TRUE)
}
names(NL21) <- NUTS1$nuts_id
qs::qsave(NL21, "data/NL21.qs")

NL21 <- rowbind(NL21, idcol = "NUTS1")
