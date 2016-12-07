### This will help you compute the shortest road trip route for 
### 1. Visiting few landmark venues in Melbourne
### 2. Visiting the top 50 cities in terms of population in Australia
### We will be using the various google APIs supported for R
### Author: Sharan Kumar Ravindran

# Setting the directory location
# Copy the data and place them in a folder named Data in the below location
setwd("C:\\Users\\skrav\\Desktop\\Deloitte\\Sharan\\Blog\\Travel\\Github")
getwd()



# package required for the fetching the geo-location based on
# name of the city
# Load the package to R
#install.packages("ggmap")
library(ggmap)

# Reading the data from source
# there are three different sets of data
# 1. venues across Victoria
# 2. Top 50 cities across Australia
locations <- read.csv("Data/TouristLocations_state.csv")
locations <- read.csv("Data/Top50Cities.csv")
head(locations)
?geocode

## Fetching the co-ordinates for the location based on the name of the place
## we use the function called geocode which is part of the package ggmap
## Create a placeholder to store the location details
## and using the loop fetch the location details for all the locations
co_ord <- geocode(as.character(locations[1,1]))
for (i in 1:nrow(locations))
{
  co_ord[i,] <- geocode(as.character(locations[i,1]))
  print(i)
}

# combine the location details with the name of the place
# that was used for searching
loc_coord <- cbind(locations$Top.Attractions , co_ord)
colnames(loc_coord) <- c("names","lon","lat")
head(loc_coord)
# saving the co-ordinates for the state
write.csv(loc_coord, "Data/Vic_places_coord.csv", row.names = FALSE)

# combine the location details with the name of the place
# that was used for searching
loc_coord <- cbind(locations$City , co_ord)
colnames(loc_coord) <- c("names","lon","lat")
head(loc_coord)
# saving the co-ordinates for top 50 cities
write.csv(loc_coord, "Data/Top50Cities_coord.csv", row.names = FALSE)

# identifying thr distance between each of them
# creating various place holders
dist <- ""
pname <- ""
p_dist <- ""
p_dist_r <- ""
p_dist_n <- ""
for (i in 1:nrow(locations))
{
  for(j in 1:nrow(locations))
  {
    pname <- cbind(as.character(locations[i,1]), as.character(locations[j,1]))
    colnames(pname) <- c("from_n", "to_n")
    dist <- mapdist(as.numeric(co_ord[i,]), as.numeric(co_ord[j,]), mode = "driving")
    p_dist <- cbind(pname, dist)
    p_dist_r <- rbind(p_dist_r, p_dist)
    print(j)
  }
  p_dist_n <- rbind(p_dist_n, p_dist_r)
  print(i)
}

# maximum of 2500 calls is allowed per day
# to check number of calls remaining use the below command
geocodeQueryCheck()

# Save the results to a file to avoid unnecessary re-run of the above code
write.csv(p_dist_n, "Data/dist_50cities.csv", row.names=FALSE)

dist_read <- read.csv("Data/dist_VIC.csv")
dist_read <- read.csv("Data/dist_50cities.csv")
# remove all the rows that has NA
dist_read_v1 <- dist_read[complete.cases(dist_read),]
head(dist_read_v1)

# choose the selected column that will be required for computing the shortest route
col_for_shrt_dist <- dist_read_v1[,c("from_n", "to_n", "km")]
head(col_for_shrt_dist)

### pivot data to suit the TSP package
### the below cast function will create a distance matrix
### where the place name will become the row and column names
### and distance between them is mentioned in the cells
library(reshape)
dist_matrix <- cast(col_for_shrt_dist, from_n ~ to_n, fun.aggregate=min)
# replace all the invalid details with zero
dist_matrix[mapply(is.infinite, dist_matrix)] <- 0
head(dist_matrix, 10)
# Convert the data into a matrix
dist_matrix <- as.matrix(dist_matrix)
dm <- dist_matrix
# compute the distance to all other place from each place
# the place with the lowest sum should be the one where we should start with
# and the place with the highest one should be avoided 
# as starting place at any cost
dm$SUM <- rowSums(dm)
startPoint <- which.min(dm$SUM) # Lowest sum
endPoint   <- which.max(dm$SUM) # Highest sum


#########TSP- to identify the ideal route
install.packages("TSP")
library(TSP)

# Converting the format of the data
atsp <- ATSP(dist_matrix)
?solve_TSP
# Compute the shourtest route
# Until the total distance is low the combination are tried
tour <- solve_TSP(atsp, method="nearest_insertion", start=startPoint)
head(tour, 10)
route <- as.data.frame(tour)
# Tour length 50 cities: TSP : 18660.06 . vic: 1368.703
## plotting the route in the order
route
names <- rownames(route)
rownames(route) <- NULL
route_name <- cbind(names, route)
route_name$order <- 1:nrow(route_name)
head(route_name)

for_direc <- merge(route_name, loc_coord, by=c("names"))
for_direc$tour <- NULL
head(for_direc)
for_direc_v1 <- for_direc[with(for_direc, order(order)), ]
for_direc_v1$order <- NULL
head(for_direc_v1)
for_direc_v1 <- rbind(for_direc_v1, for_direc_v1[1,])
write.csv(for_direc_v1, "Data\\short_route_vic_tsp.csv", row.names = FALSE)
# Plotting the cities
# Leaflet is an package used for making the plot on the google maps
install.packages("leaflet")
library(leaflet)

# Ploting the top 50 cities in the map
leaflet() %>%
  addTiles() %>%  
  setView(134.118785, -27.289612, zoom = 4) %>%
  addMarkers(data = for_direc_v1, lat = ~lat, lng = ~lon)

# Ploting the venue across victoria in the map
leaflet() %>%
  addTiles() %>%  
  setView(134.118785, -27.289612, zoom = 4) %>%
  addMarkers(data = for_direc_v1, lat = ~lat, lng = ~lon)


# plotting the route for Australia trip

leaflet() %>%
  addTiles() %>%
  addPolylines(data = for_direc_v1, lat = ~lat, lng = ~lon)

########################################################
##########################Using Genetic Algorithm
#######################################################
install.packages("GA")
library(GA)

# Creating a distance matrix suitable for the GA Algorithm
dist_matrix <- as.data.frame(dist_matrix)
rownames(dist_matrix) <- dist_matrix$from_n
dist_matrix$from_n <- NULL
dist <- dist_matrix
D <- as.matrix(dist)


#Function to calculate tour length 

tourLength <- function(tour, distMatrix) {
  tour <- c(tour, tour[1])
  route <- embed(tour, 2)[,2:1]
  sum(distMatrix[route])
}

#Firness function to be maximized

tspFitness <- function(tour, ...) 1/tourLength(tour, ...)

####### model iterations for VIC trip
GA <- ga(type = "permutation", fitness = tspFitness, distMatrix = D,
         min = 1, max = 10, popSize = 10, maxiter = 40000,
         run = 20000, pmutation = 0.3)

####### model iterations for Australia trip
GA <- ga(type = "permutation", fitness = tspFitness, distMatrix = D,
         min = 1, max = 50, popSize = 50, maxiter = 500000,
         run = 50000, pmutation = 0.4)

head(D)
# Model details
summary(GA)
tour_ga <- GA@solution[1, ]
tour_ga <- c(tour_ga, tour_ga[1])
route_ga <- embed(tour_ga, 2)[,2:1]
# Distance of shortest route
sum(dist[route_ga])
# Distance by GA- 50 cities: 18123.69 . vic: 1323.281

# Optimal route with matrix position
route_ga <- as.data.frame(route_ga)
colnames(route_ga) <- c("from", "to")
route_ga$order <- 1:nrow(route_ga)
head(route_ga)

# Replacing with the name of the places
destination <- rownames(dist_matrix)
destination <- as.data.frame(destination)
destination$from <- 1:nrow(destination)
head(destination)
for_direc_ga_v2 <- merge(route_ga, destination, by=c("from"))
for_direc_ga_v2_1 <- for_direc_ga_v2[with(for_direc_ga_v2, order(order)), ]
head(for_direc_ga_v2_1)

# To plot the actual route
path <- data.frame(for_direc_ga_v2_1$destination, for_direc_ga_v2_1$order)
colnames(path) <- c("names", "order")
path_full <- merge(path, loc_coord, by=c("names"))
path_full <- path_full[with(path_full, order(order)), ]
path_full$order <- NULL
path_full <- rbind(path_full, path_full[1,])
head(path_full)
write.csv(path_full, "Data\\short_route_aus_ga_v1.csv", row.names = FALSE)
# plotting the route

library(leaflet)
leaflet() %>%
  addTiles() %>%
  addPolylines(data = path_full, lat = ~lat, lng = ~lon)

