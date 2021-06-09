#Yearly Average Runoff

library(ggplot2)
library(tidyr)
library(tidyverse)
install.packages("gifski")
library(gganimate)
library(gifski)
library(hrbrthemes)

#nc files from https://disc.gsfc.nasa.gov/datasets/FLDAS_NOAH01_C_GL_M_001/summary?keywords=FLDAS and in Box

#extract_nc() function defined in functions_created.R script used to obtain list of values for interested variable
#list_to_df() function defined in functions_creaed.R script used to obtain dataframe of values for interested variable from list

files <- paste("ges_disc/FLDAS_NOAH01_C_GL_M.A", rep(1982:2020, each = 12), 
               c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"), ".001.nc.SUB.nc4", sep = "")
files <- c(files, paste("ges_disc/FLDAS_NOAH01_C_GL_M.A", rep(2021, 3), 
                        c("01", "02"), ".001.nc.SUB.nc4", sep = ""))


shp <- st_read("tana_outline.shp") #Shapefile of the Upper Tana River Basin

## Surface runoff Animated
surfrun_list <- extract_nc(files, "Qs_tavg", start = c(start_lon, start_lat, 1), count = c(count_lon, count_lat, 1))
surfrun_df <- list_to_df(surfrun_list)

surfrun_df_group <- group_by(surfrun_df, Year, Longitude, Latitude)
avgyear_surfrun <- summarise(surfrun_df_group, "Surface Runoff (kg m-2 s-1)" = mean(Vals))
ungroup(surfrun_df_group)

surf_run_p <- ggplot() +  
  geom_tile(avgyear_surfrun, mapping = aes(x = Longitude, y = Latitude, fill = `Surface Runoff (kg m-2 s-1)`)) +
  scale_fill_gradient2(limits = c(1.410458e-07, 4.42616e-05)) +
  geom_sf(data=shp,fill = "transparent", col = "black") 
surf_run_p2 <- surf_run_p + 
  labs(title = "Average Surface Runoff ({frame_time})") +
  transition_time(avgyear_surfrun$Year)
animate(surf_run_p2, nframes = 150, end_pause = 5, renderer = gifski_renderer("surfrunoff_yearly.gif"))

## Surfrunn198590-201520
surfrun_df_8590 <- filter(surfrun_df, Year >= 1985 & Year <= 1990)
surfrun_8590group <- group_by(surfrun_df_8590, Longitude, Latitude)
surfrun8590_avg <- summarise(surfrun_8590group, "1985-1990 Average Surface Runoff (mm/day)" = mean(Vals) * 86400)
ungroup(surfrun_8590group)

ggplot() + 
  geom_tile(surfrun8590_avg, mapping = aes(x = Longitude, y = Latitude, fill = `1985-1990 Average Surface Runoff (mm/day)`)) +
  scale_fill_gradient2(limits = c(0, 2.5)) +
  geom_sf(data = shp, fill = "transparent", col = "black")

## Precipitation Animated
precip_list <- extract_nc(files, "Rainf_f_tavg", start = c(start_lon, start_lat, 1), count = c(count_lon, count_lat, 1))
precip_df <- list_to_df(precip_list)

precip_df_group <- group_by(precip_df, Year, Longitude, Latitude)
avgyear_precip <- summarise(precip_df_group, "Precipitation Rate (mm/day)" = mean(Vals) * 86400)
ungroup(precip_df_group)

precip_p <- ggplot() +
  geom_tile(avgyear_precip, mapping = aes(x = Longitude, y = Latitude, fill = `Precipitation Rate (mm/day)`)) +
  scale_fill_gradient2(low = "white", mid = "skyblue", high = "midnightblue", limits = c(0, 15), midpoint = 8) +
  geom_sf(data=shp,fill = "transparent", col = "black")
precip_p2 <- precip_p +
  labs(title = "Average Precipitation Rate ({frame_time})") +
  transition_time(avgyear_precip$Year)
animate(precip_p2, nframes = 150, end_pause = 5, renderer = gifski_renderer("preciprate_yearly.gif"))


## Precipitation Rainy Season Animated
precip_rainy_df <- filter(precip_df, Month == 1 | Month <= 5 & Month >= 3 | Month <= 12 & Month >= 10)

precip_rainy_group <- group_by(precip_rainy_df, Year, Longitude, Latitude)
avgyr_precip_rainy <- summarise(precip_rainy_group, "Precipitation Rate (mm/day)" = mean(Vals) * 86400)
ungroup(precip_rainy_group)


precip_rain_p <- ggplot() +
  geom_tile(avgyr_precip_rainy, mapping = aes(x = Longitude, y = Latitude, fill = `Precipitation Rate (mm/day)`)) +
  scale_fill_gradient2(low = "white", mid = "skyblue", high = "midnightblue", limits = c(0, 15), midpoint = 8) +
  geom_sf(data=shp,fill = "transparent", col = "black")
precip_rain_p2 <- precip_rain_p +
  labs(title = "Average Rainy Season Precipitation Rate ({frame_time})") +
  transition_time(avgyr_precip_rainy$Year)
animate(precip_rain_p2, nframes = 150, end_pause = 5, renderer = gifski_renderer("preciprate_rainy_yearly.gif"))

## Precipitation Dry Season Animated
precip_dry_df <- filter(precip_df, Month == 2 | Month <= 9 & Month >= 6)

precip_dry_group <- group_by(precip_dry_df, Year, Longitude, Latitude)
avgyr_precip_dry <- summarise(precip_dry_group, "Precipitation Rate (mm/day)" = mean(Vals) * 86400)
ungroup(precip_dry_group)

precip_dry_p <- ggplot() +
  geom_tile(avgyr_precip_dry, mapping = aes(x = Longitude, y = Latitude, fill = `Precipitation Rate (mm/day)`)) +
  scale_fill_gradient2(low = "white", mid = "skyblue", high = "midnightblue", limits = c(0, 15), midpoint = 8) +
  geom_sf(data=shp,fill = "transparent", col = "black")
precip_dry_p2 <- precip_dry_p +
  labs(title = "Average Dry Season Precipitation Rate ({frame_time})") +
  transition_time(avgyr_precip_dry$Year)
animate(precip_dry_p2, nframes = 150, end_pause = 5, renderer = gifski_renderer("preciprate_dry_yearly.gif"))


## Soil Moisture
smc_list <- extract_nc(files, "SoilMoi00_10cm_tavg", start = c(start_lon, start_lat, 1), count = c(count_lon, count_lat, 1)) 
smc_df <- list_to_df(smc_list)

smc_group <- group_by(smc_df, Year, Longitude, Latitude)
avgyr_smc <- summarise(smc_group, "Soil Moisture (m3/m3)" = mean(Vals))
ungroup(smc_group)

smc_p <- ggplot() +
  geom_tile(avgyr_smc, mapping = aes(x = Longitude, y = Latitude, fill = `Soil Moisture (m3/m3)`)) +
  scale_fill_gradient2(low = "white", mid = "lightsteelblue", high = "midnightblue", limits = c(0.1, 0.42), midpoint = 0.25) +
  geom_sf(data=shp,fill = "transparent", col = "black")
smc_p2 <- smc_p +
  labs(title = "Average Soil Moisture Content at 0-10cm ({frame_time})") +
  transition_time(avgyr_smc$Year)
animate(smc_p2, nframes = 150, end_pause = 5, renderer = gifski_renderer("smc0010_yearly.gif"))


