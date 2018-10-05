# load libraries
require(tidyverse)
require(sf)
source("R/mpmm.r")
source("R/mapr.r")

# import data
d <- read.csv("dat/hp2_ctrw.csv", stringsAsFactors = FALSE) %>%
  tbl_df() %>%
  select(-X) %>%
  rename_all(tolower) %>%
#  filter(id == "hp5-L764-18") %>%
  mutate(x = x/1000, y = y/1000) %>%
  filter(av.dur > 0) %>%
  mutate(av.depth = ifelse(av.depth >=8 , av.depth, NA)) %>%
  mutate(date = lubridate::ymd_hms(date, tz = "UTC")) %>%
  filter(id != "hp2-9320-04", id != "hp2-9349-04")

# plot data to check and demonstrate new 'mapr' function
d_sf <- st_as_sf(d, coords = c("lon", "lat")) %>%
  st_set_crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

prj = "+proj=laea +lat_0=75 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
world_shp <- mapr(d, prj, buff = 500000)

ggplot() +
  geom_sf(aes(), data = world_shp) + 
  geom_sf(aes(), data = d_sf)

# run lambda mpmm model
fit <- mpmm(~ 1 + (1 | id), d, verbose = TRUE)

# check convergence
fit$opt

# plot using sf
fit_sf <- st_as_sf(fit$data, coords = c("lon", "lat")) %>%
  st_set_crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

# plot output
ggplot() + 
  geom_sf(aes(colour = fit$fitted$gamma), size = 0.5, data = fit_sf) +
  scale_colour_viridis_c() +
  facet_wrap(~id)

ggplot() + 
  geom_sf(aes(colour = fit$fitted$lambda), size = 0.5, data = fit_sf) +
  scale_colour_viridis_c() +
  facet_wrap(~id)

ggplot() + 
  geom_path(aes(x = lambda, y = gamma, group = id), data = fit$fitted) +
  facet_wrap(~id)
