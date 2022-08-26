library(sf)
library(dplyr)

nc  = st_read("https://www2.stat.duke.edu/~cr173/Sta444_Fa18/slides/data/nc_counties.gpkg", quiet=TRUE, stringsAsFactors=FALSE)
air = st_read("https://www2.stat.duke.edu/~cr173/Sta444_Fa18/slides/data/airports.gpkg", quiet=TRUE, stringsAsFactors=FALSE)
hwy = st_read("https://www2.stat.duke.edu/~cr173/Sta444_Fa18/slides/data/us_interstates.gpkg", quiet=TRUE, stringsAsFactors=FALSE)

## What counties are adjacent to Durham County

neighbors = st_touches(
  filter(nc, COUNTY == "Durham County"),
  nc
) %>% 
  unlist()

plot(st_geometry(nc))
plot(st_geometry(nc[neighbors,]), add=TRUE, col="lightblue")
plot(st_geometry(nc[nc$COUNTY == "Durham County",]), add=TRUE, col="red")



## What counties have more than 4 neighbors

library(purrr)

many_neigh = st_touches(nc) %>% map_int(length) %>% {which(. > 4)} 

plot(st_geometry(nc))
plot(st_geometry(nc[many_neigh,]), add=TRUE, col="lightblue")

## What counties have an airport

nc_air = st_intersection(nc, air)

have_airport = st_intersects(nc, nc_air) %>% map_int(length) %>% {which(. > 0)}

plot(st_geometry(nc))
plot(st_geometry(nc[have_airport,]), add=TRUE, col="lightblue")
plot(nc_air, add=TRUE, col="red", pch=16)
