### Ant colony optimisation for TSP

# Import TSP library for TSP benchmarks, NN and 2-opt algorithms

library(TSP)
library(ggplot2)

# Calculate attractness and return probability distribution

attract <- function(pheromone.drop){
  n_i <- 1/dist[,]
  n_i[n_i == Inf] <- 0
  prob <- (n_i^beta)*(pheromone.drop[,]^alpha)/(sum((n_i^beta)*(pheromone.drop[,]^alpha))) # Compute probabilities
  return (prob)
}


# Calculate route for one ant 

route <- function(start, pheromone.drop, theta){
  # Parameters may follow a normal distribution instead of only one value
  #alpha <- rnorm(1, 1, 0.1)
  #beta <- rnorm(1, 3, 0.5)
  idx <- 1:n
  idx <- idx[-start]
  i <- start
  route <- data.frame()
  route <- rbind(route, start)
  shannon <- 0
  for (k in 1:(n-1)){
    if (length(idx) == 1){
      route <- rbind(route, idx)
    } else {
      
      if (isSha){
        shannon <- shannon - sum(log(theta[idx, i])*theta[idx, i])
      }
      
      choice <- sample(idx, 1, prob = theta[idx, i]) # Choose one next node
      idx <- idx[-which(idx == choice)] # Delete this node for next path
      route <- rbind(route, choice) # Adding the node to the route
      i <- choice # Assign new starting node
    }
  }
  
  route <- rbind(route, start) # Make it a cycle
  if (isSha){
    shannon <- shannon/(n-1)
    sha[length(sha)+1] <<- shannon 
  }

  return (route)
}

# Calculate route cost

tour_cost <- function (tour){
  sum <- 0
  for (k in 1:(n)){
    i <- tour[k,1]
    j <- tour[k+1,1]
    sum <- sum + dist[i,j] 
  }
  return (sum)
}

# Return route and cost for each ants

launch <- function (elit.tour, pheromone){
  ants_tour.cost <- data.frame()
  ants_tour.route <- matrix()
  pheromone.drop <- pheromone
  theta <- attract(pheromone.drop)
  for (k in 1:m){
    
    if (isDrop){
      pheromone.drop <- drop_out(elit.tour, pheromone)
    }
    
    tour <- route(ceiling(runif(1, 0, n)), pheromone.drop, theta)
    cost <- tour_cost(tour)
    ants_tour.cost <- rbind(ants_tour.cost, cost)
    ants_tour.route <- cbind(ants_tour.route, tour)
  }
  return (c(ants_tour.route[-1], ants_tour.cost))
}

# Update pheromone matrix (basic)

update <- function (ants_tour.route, ants_tour.cost, pheromone){
  pheromone[,] <-  (1-evaporation)*pheromone[,]
  for (k in 1:m){
    for (l in 1:(n)){
      i <- ants_tour.route[l,k]
      j <- ants_tour.route[l+1,k]
      pheromone[i,j] <- pheromone[i,j] + 1/ants_tour.cost[k,1]
    }
  }
  return (pheromone)
}

# Update pheromone matrix (elitist)

update_elitist <- function (elit, ants_tour.route, ants_tour.cost, i, pheromone){
  pheromone[,] <-  (1-evaporation)*pheromone[,]
    for (l in 1:(n)){
      i <- ants_tour.route[l,elit]
      j <- ants_tour.route[l+1,elit]

      # Drop Connect
      
      prob <- 0.1
      drop <- sample(0:1,1, prob = c(prob,1-prob))
      if (isDrop) {
        pheromone[i,j] <- pheromone[i,j] + e/ants_tour.cost[elit,1]*drop
        #pheromone[i,j] <- pheromone[i,j] + e/ants_tour.cost[elit,1]*rnorm(1,1,0.2)
      } else {
        pheromone[i,j] <- pheromone[i,j] + e/ants_tour.cost[elit,1]
      }
      
    }
    for (k in 1:m){
      for (l in 1:(n)){
        i <- ants_tour.route[l,k]
        j <- ants_tour.route[l+1,k]
        pheromone[i,j] <- pheromone[i,j] + 1/ants_tour.cost[k,1]
    }
    }
  return (pheromone)
}


# Drop Out Function (randomly ban one interesting arc)

drop_out <- function(tour, pheromone){
    choice <- sample(1:n,1)
    i <- tour[choice]
    j <- tour[choice+1]
    pheromone.drop <- pheromone
    pheromone.drop[i, j] <- 0
    return (pheromone.drop)
}


# Launch ACO function

aco <- function (){
  means <- matrix(ncol = 1, nrow = k)
  min <- Inf
  #best.tour <- c(1:n, 1)
  elit.tour <- c(1:n, 1)
  pheromone <- matrix(nrow = n, ncol = n)
  pheromone[,] <- t_0 # Inititialize pheromone 
  for (i in 1:k){
    results <- launch(elit.tour, pheromone)
    ants_tour.route <- as.data.frame(results[1:m])
    ants_tour.cost <- as.data.frame(results[[m+1]])
    elit <- which.min(ants_tour.cost[,1]) # Choose elit for ASE
    elit.tour <- ants_tour.route[,elit]
    
    if (min > ants_tour.cost[elit,1]){
      min <- ants_tour.cost[elit,1]
      best.elit <- elit
      best.tour <- elit.tour
    }
    
    if (isElit){
      pheromone <- update_elitist(elit, ants_tour.route, ants_tour.cost, i, pheromone) # Update pheromone for ASE
    } else {
      pheromone <- update(ants_tour.route, ants_tour.cost, pheromone) # Update pheromone for AS
    }
    
    #if ((i%%100) == 0){
      # Open a pdf file
      #pdf(paste("rplot",i,".pdf")) 
      #heatmap.2(pheromone,dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none')
      # Close the pdf file
      #dev.off() 
    #}

    #cost.sd[i] <<- sd(ants_tour.cost[,1])/mean(ants_tour.cost[,1])
    means[i] <-  min(ants_tour.cost[,1])
  }
  return (c(means, best.tour))
}



sha <- matrix()
tsp <- TSP(dist)
best_nn <- solve_TSP(tsp, method = "repetitive_nn", two_opt = FALSE)
opt <- tour_length(best_nn)
n <- nrow(dist) # Number of cities
m <- n # Number of ants (same as number of cities)
cost.sd <-  matrix(nrow = n, ncol = 1)
C.nn <- opt # Nearest neighbor route
e <- n # For elitist system
evaporation <- 0.5 # Coefficient of evaporation
alpha <- 1
beta <- 3.5
k <- 100 # Number of iterations
isElit <- TRUE
isDrop <- TRUE
isSha <- FALSE

if (isElit){
  t_0 <- (e+m)/(evaporation*C.nn) # Initial pheromone for ASE
} else {
  t_0 <- m/C.nn # Initial pheromone value for AS
}

test <- aco()


## Analysis part (not useful if you're only interesting in the algorithm code)

# Some variables useful for my analysis

plot(1:k, test[1:k], xlab = "Number of iterations", ylab= "Tour length", ylim = c(110000, 150000))
plot(1:k, cost.sd[1:k], xlab = "Number of iterations", ylab= "Coefficient of variation", ylim = c(0.04, 0.10))
test.as.1000 <- test
tour.as.1000 <- test[1001:1051]
cv.as.1000 <- cost.sd
test.ase.1000 <- test
tour.ase.1000 <- test[1001:1051]
cv.ase.1000 <- cost.sd
test.ased.1000 <- test
tour.ased.1000 <- test[1001:1051]
cv.ased.1000 <- cost.sd


# Calculate shannon index

sha.mean <- matrix(ncol=1, nrow = k)
for (i in 0:(k-1)){
  sha.mean[i+1] <- mean(sha[(i*50+2):(52+i*50)])
}
plot(1:(k-1), sha.mean[1:(k-1)], xlab = "Number of iterations", ylab = "Shannon index", ylim = c(0.022,0.026))
sha.as.1000 <- sha.mean
sha.ase.1000 <- sha.mean
sha.ased.1000 <- sha.mean

# Run several experiments

n.exp <- 10 # Number of experiments
aco.10 <- replicate(n.exp, aco())

# 

df <- as.data.frame(aco.10)
opt <- 675

best.first <- 1:n.exp
for (i in 1:n.exp){
  best.first[i] <- which(ceiling(df[,i]) == 7545)[1]
}

best <- 1:n.exp
for (i in 1:n.exp){
  best[i] <- min(df[,i])
}

df.sd <- 1:k
for (i in 1:k){
  df.sd[i] <- sd(df[i,])
}

x <- 1:k
y <- rowSums(df)/n.exp

qplot(x,y)+geom_errorbar(aes(x=x, ymin=y-df.sd, ymax=y+df.sd), width=0.5) + geom_hline(yintercept = opt, colour = "red")

min(best)
max(best)
mean(best)
sd(best)
error <- ((mean(best)-opt)/opt)*100
error


# Load TSP data

tsp <- read_TSPLIB("berlin52.tsp")
tsp <- read_TSPLIB("eil76.tsp") # opt = 538
tsp <- read_TSPLIB("kroB200.tsp")
tsp <- read_TSPLIB("eil51.tsp") #opt = 426
tsp <- read_TSPLIB("st70.tsp") # opt = 675
matrix <- as.matrix(tsp)
dist <- dist(tsp, method = "euclidian")
dist <- as.matrix(dist_matrix)

test <- solve_TSP(tsp, method = "two_opt") # 8182
test <- solve_TSP(tsp, method = "two_opt", tour = 1:52)
tour_length(test)
tour <- as.matrix(test)


cost <- matrix(ncol = 1, nrow = 10)
for (i in 1:10){
  test <- solve_TSP(tsp, method = "nn", two_opt = TRUE)
  cost[i] <- tour_length(test)
}
min(cost)
max(cost)
mean(cost)
sd(cost)

