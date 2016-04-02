library(foreach)
library(GenSA)

tsp_sa <- function(max_iter, in_iter, points, temp,
                   lowest_temp=0, consecutive=false)# Solve tsp problem with simulated annealing
                                                    # temp is initial temperature
                                                    # Points is data_frame generated with generate_data
                                                    # in_iter is number of iterations for fixed temperature
{
  costlist <- rep(0, max_iter) # Keep track of cost value in each iteration
  templist <- rep(0, max_iter) # Keep track of temperature changes
  temp_factor <- (temp/lowest_temp)^(-1/max_iter)
  numaccept <- 0  # How many points were swapped
  
  for(i in 1:max_iter) {
    costlist[i] <- distance(points) # Save current cost
    templist[i] <- temp # Save current temperature
    
    for(t in 1:in_iter) {
      
      if(consecutive) { # Consecutive swap
        a = sample(1:nrow(points), 1) # Randomly choose one point for swapping with his successor
        
        if(a==nrow(points)) {
          a <- c(a, 1)
        }
        else {
          a <- c(a,a+1)
        }
      }
      else { # Arbitrary swap
        a = sample(1:nrow(points), 2) # Randomly choose two points for swapping
      }
      tmp = points # Create copy of points to compare distance after swapping
      tmp[a,] = tmp[rev(a),] # Swap points
      
      U = runif(1)
      if(U < exp((distance(points)-distance(tmp))/temp)) { # Decide whether to keep swapped points
                                                           # based on temp and random number U
        points <- tmp
        numaccept <- numaccept + 1 
      }
    }
    temp <- temp * 0.999 # Update temperature
    temp <-temp*temp_factor 
  }
  return(list(points, costlist, templist, numaccept))
}
generate_data <- function(size, seed=1, max=1) # Generate random set of 2D points
{
  set.seed(seed)
  x <- runif(size)
  y <- runif(size)
  
  df <- data.frame(x,y)
  return (df*max)
}

euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2)) # Euclidean distance between
                                                      #two points

distance <- function(X) # Calculates distance between all
                        #neighbouring points ((i,i+1) and (first,last) ) in X 
{


  sum(foreach(i = 1:(nrow(X)-1), .combine = c ) %do% euc.dist(X[i,],X[i+1,]),
             euc.dist(X[nrow(X),],X[1,]))

}

draw_path <- function(X) # Draws how the points are visited
{
  plot(X)
  for(i in 1:(nrow(X)-1)) {
    segments(X[i,1],X[i,2], X[i+1,1], X[i+1,2]) # Draw line between consecutive points
  }
  segments(X[1,1],X[1,2], X[nrow(X),1], X[nrow(X),2]) # Draw line between first and last
}

draw_changes <- function(changelist) # Draws the changes of temperature or distance
{
  plot(seq(1:length(changelist)), changelist,type='o' , pch=10, cex=.2)
}