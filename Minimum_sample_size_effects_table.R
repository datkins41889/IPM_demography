data_out = sf::st_drop_geometry(data_out)
View(data_out)

library(tidyverse)

data_out$Site = gsub("Ft Valley - COC-S1A", "FVEF S1", data_out$Site)
data_out$Site = gsub("Ft Valley - COC-S1B", "FVEF S1", data_out$Site)
data_out$Site = gsub("Ft Valley - COC-S2A", "FVEF S2", data_out$Site)
data_out$Site = gsub("Ft Valley - COC-S2B", "FVEF S2", data_out$Site)
data_out$Site = gsub("Ft Valley - COC-S3A", "FVEF S3", data_out$Site)

test = data_out %>%
  group_by(species, Site, z_Year) %>%
  summarise(n = n())
test


nrow(test)
Nzero = test %>%
  group_by(species, Site) %>%
  summarise(n=n())
nrow(Nzero)
Nzero %>%
  group_by(species) %>%
  summarise(populations = length(n), pop_years = sum(n))

nrow(test[test$n >= 5,])
Nfive = test[test$n >= 5,] %>%
  group_by(species, Site) %>%
  summarise(n=n())
nrow(Nfive)
Nfive %>%
  group_by(species) %>%
  summarise(populations = length(n), pop_years = sum(n))

nrow(test[test$n >= 10,])
Nten = test[test$n >= 10,] %>%
  group_by(species, Site) %>%
  summarise(n=n())
nrow(Nten)
Nten %>%
  group_by(species) %>%
  summarise(populations = length(n), pop_years = sum(n))

nrow(test[test$n >= 15,])
Nfifteen = test[test$n >= 15,] %>%
  group_by(species, Site) %>%
  summarise(n=n())
nrow(Nfifteen)
Nfifteen %>%
  group_by(species) %>%
  summarise(populations = length(n), pop_years = sum(n))

nrow(test[test$n >= 20,])
Ntwenty = test[test$n >= 20,] %>%
  group_by(species, Site) %>%
  summarise(n=n())
nrow(Ntwenty)
Ntwenty %>%
  group_by(species) %>%
  summarise(populations = length(n), pop_years = sum(n))

nrow(test[test$n >= 30,])
Nthirty = test[test$n >= 30,] %>%
  group_by(species, Site) %>%
  summarise(n=n())
nrow(Nthirty)
Nthirty %>%
  group_by(species) %>%
  summarise(populations = length(n), pop_years = sum(n))




