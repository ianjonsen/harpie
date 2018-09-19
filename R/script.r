require(tidyverse)

d <- read.csv("dat/hp2_ctrw.csv", stringsAsFactors = FALSE) %>%
  tbl_df() %>%
  select(-X) %>%
  rename_all(tolower) %>%
#  filter(id == "hp5-L764-18") %>%
  mutate(x = x/1000, y = y/1000) %>%
  filter(av.dur > 0) %>%
  mutate(av.depth = ifelse(av.depth >=8 , av.depth, NA)) %>%
  mutate(date = lubridate::ymd_hms(date, tz = "UTC")) %>%
  filter(id != "hp2-9349-04", id != "hp2-9320-04")


fit <- mpmm(~ 1 + (1 | id), d, verbose = TRUE)
