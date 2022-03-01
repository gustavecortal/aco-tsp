### Prepare data and plot the best path on Google Maps

## Import libraries

library(dplyr)
library(purrr)
library(ggmap)
library(gmapsdistance)
library(TSP)
library(ggplot2)

## Prepare data

# Google API

api_key <- "AIzaSyC4aUz9VbYr3NoiMfIka_UAEpp4Y8rJW9w"
set.api.key(api_key)
register_google(api_key)

# Import data

data <- read.csv(file = "ubike_tp.csv", header = TRUE)

# Select random stations and keep latitude and longitude values

cities.idx <- sample(400, 50)
data <- data[cities.idx,]
latlon <- paste(data$latitude,data$longitude, sep = '+')

# Compute distance matrix using Google API 

distances <- gmapsdistance(origin = latlon[1:50],
                           destination = latlon[1:50],
                           combinations = "all", mode = "driving")

# Convert to Distance matrix

dist_matrix <- as.dist(distances$Distance[,-1])

# Apply NN+2-opt

tsp <- TSP(dist_matrix)
tour <- solve_TSP(tsp, "repetitive_nn", two_opt = TRUE)
tour_order <- as.integer(tour)
tour_length(tour)

# Return route order

new_laton_order <- latlon[tour_order]
new_laton_order[51] <- new_laton_order[1]
cities.idx[tour_order] 

## Compute route order in Google Maps and plot the map

map <- get_map(location = c(lon = 121.55, lat = 25.070), zoom = 12, maptype = "roadmap")
map <- get_map(location = "taipei", maptype = "roadmap", zoom = 11)

from <- new_laton_order[1]
to <- new_laton_order[2]
route_df <- route(from, to, structure = "route")

for(i in 2:50){
  from <- new_laton_order[i]
  to <- new_laton_order[i+1]
  route_df <- rbind(route_df, route(from, to, structure = "route"))
}

ggmap(map, extent = "device") +
  geom_point(data = data, aes(x = longitude, y = latitude),
             color = 'blue', size = 3.5) +
  geom_path(
    aes(x = lon, y = lat),  colour = "red", size = 1.2,
    data = route_df, lineend = "round"
  )

# Some interestings results on the selected path

sum(as.matrix(route_df[2]), na.rm = TRUE) # Tour length (m)
sum(as.matrix(route_df[3]), na.rm = TRUE) # Tour length (km)
sum(as.matrix(route_df[6]), na.rm = TRUE) # Minutes
