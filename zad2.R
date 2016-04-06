library(stats)
library(ConsPlan)

generate_img <- function(N, delta=0.5) # Generates an NxN binary image
                                       # with density of black points equal to delta
{
  m <- rep(0, N^2) # Vector of zeros
  x <- sample(1:(N^2), delta*N^2) # Pick random values to blacken
  m[x] <- 1 # Blacken chosen values
  
  m <- matrix(m, ncol=N) # Create matrix from vector
  image(m, axes = FALSE, col=grey(seq(1,0))) # Display generated image
  return(m)
}


allsame_cost <- function(X) # Returns the sum of numbers of neighbours
                            # with different color for each pixel
{
  X <- matrix(X, ncol=sqrt(length(X))) # sann passes X as list
  sum <- 0
  for(i in 2:(nrow(X)-1)) { # Regular case - 8 neighbours
    for(j in 2:(ncol(X)-1)) {
      tmp <- sum(X[i-1,j-1], X[i,j-1], X[i+1,j-1], X[i-1,j], X[i+1,j], 
                 X[i-1,j+1], X[i,j+1], X[i+1,j+1])
      sum <- sum + abs(tmp-8*X[i,j])
    }
  }
  # 4 Corner cases - 3 neighbours
  
  sum <- sum + abs(3*X[1,1]-sum(X[1,2],X[2,1],X[2,2])) # Top-left
  sum <- sum + abs(3*X[nrow(X),1]-sum(X[nrow(X),2],X[nrow(X)-1,1],X[nrow(X)-1,2])) # Bottom-left
  sum <- sum + abs(3*X[nrow(X),ncol(X)]-sum(X[nrow(X),ncol(X)-1],X[nrow(X)-1,ncol(X)],X[nrow(X)-1, ncol(X)-1])) #Bottom-right
  sum <- sum + abs(3*X[1,ncol(X)]-sum(X[1,ncol(X)-1],X[2,ncol(X)],X[2,ncol(X)-1])) # Top-right
  
  # 4 Side cases - 5 neighbours
  
  for(i in 2:(ncol(X)-1)) {
    sum <- sum + abs(5*X[1,i] - sum(X[1,i-1], X[2,i-1], X[2,i], X[2,i+1], X[1,i+1])) # Top
    sum <- sum + abs(5*X[i,1] - sum(X[i-1,1], X[i-1,2], X[i,2], X[i+1,2], X[i+1,1])) # Left
    sum <- sum + abs(5*X[i,ncol(X)] - sum(X[i-1, ncol(X)], X[i-1, ncol(X)-1], X[i, ncol(X)-1], X[i+1, ncol(X)-1], X[i+1, ncol(X)])) # Right
    sum <- sum + abs(5*X[nrow(X),i] - sum(X[nrow(X),i-1], X[nrow(X)-1,i-1], X[nrow(X)-1,i], X[nrow(X)-1,i+1], X[nrow(X),i+1])) # Bottom
  }
  return (sum)
}

allsame_next <- function(X) # Changes bottom-left neighbours to white and top-right to black
{
  X <- matrix(X, sqrt(length(X)))
  
  x <- sample(2:(nrow(X)-1), 1)
  y <- sample(2:(ncol(X)-1), 1)
  
  val1 <- 0
  val2 <- 1

  X[x-1,y] <- val1
  X[x-1,y-1] <- val1
  X[x-1,y+1] <- val1
  X[x,y-1] <- val1
  
  X[x,y+1] <- val2
  X[x+1,y] <- val2
  X[x+1,y-1] <- val2
  X[x+1,y+1] <- val2
  
  return(c(X))
}

equality_cost <- function(X) # Cost function returning the sum of the absolute differences
                             # between number of black and white neighbours for each pixel
                             # E.g. pixel may have up to 8 neighbours:
                             # if it has 3 black neighbours and 5 white neighbours
                             # then it returns |3-5|=2
{
  X <- matrix(X, ncol=sqrt(length(X))) # sann passes X as list
  
  sum <- 0
  for(i in 2:(nrow(X)-1)) { # Regular case - 8 neighbours
    for(j in 2:(ncol(X)-1)) {
      tmp <- sum(X[i-1,j-1], X[i,j-1], X[i+1,j-1], X[i-1,j], X[i+1,j], # if tmp is 4
                       X[i-1,j+1], X[i,j+1], X[i+1,j+1])               # then cost should be 0
      sum <- sum + abs(tmp-4)
      
    }
  }
  # 4 Corner cases - 3 neighbours
  
  sum <- sum + abs(1.5-sum(X[1,2],X[2,1],X[2,2])) # Top-left
  sum <- sum + abs(1.5-sum(X[nrow(X),2],X[nrow(X)-1,1],X[nrow(X)-1,2])) # Bottom-left
  sum <- sum + abs(1.5-sum(X[nrow(X),ncol(X)-1],X[nrow(X)-1,ncol(X)],X[nrow(X)-1, ncol(X)-1])) #Bottom-right
  sum <- sum + abs(1.5-sum(X[1,ncol(X)-1],X[2,ncol(X)],X[2,ncol(X)-1])) # Top-right
  
  # 4 Side cases - 5 neighbours
  
  for(i in 2:(ncol(X)-1)) {
    sum <- sum + abs(2.5 - sum(X[1,i-1], X[2,i-1], X[2,i], X[2,i+1], X[1,i+1])) # Top
    sum <- sum + abs(2.5 - sum(X[i-1,1], X[i-1,2], X[i,2], X[i+1,2], X[i+1,1])) # Left
    sum <- sum + abs(2.5 - sum(X[i-1, ncol(X)], X[i-1, ncol(X)-1], X[i, ncol(X)-1], X[i+1, ncol(X)-1], X[i+1, ncol(X)])) # Right
    sum <- sum + abs(2.5 - sum(X[nrow(X),i-1], X[nrow(X)-1,i-1], X[nrow(X)-1,i], X[nrow(X)-1,i+1], X[nrow(X),i+1])) # Bottom
  }
  return (sum)
  
}

equality_next <- function(X) # Function generating next state 
                             # Toggles the value of each neighbour of a randomly chosen element
                             # (except for border elements)
{
  X <- matrix(X, sqrt(length(X)))
  
  x <- sample(2:(nrow(X)-1), 1)
  y <- sample(2:(ncol(X)-1), 1)
  
  if(X[x-1,y-1]==1) {
    X[x-1,y-1]=0
  } else {
    X[x-1,y-1]=1
  }
  
  if(X[x-1,y]==1) {
    X[x-1,y]=0
  } else {
    X[x-1,y]=1
  }
  
  if(X[x-1,y+1]==1) {
    X[x-1,y+1]=0
  } else {
    X[x-1,y+1]=1
  }
  
  if(X[x,y-1]==1) {
    X[x,y-1]=0
  } else {
    X[x,y-1]=1
  }
  
  if(X[x,y+1]==1) {
    X[x,y+1]=0
  } else {
    X[x,y+1]=1
  }
  
  if(X[x+1,y-1]==1) {
    X[x+1,y-1]=0
  } else {
    X[x+1,y-1]=1
  }
  
  if(X[x+1,y]==1) {
    X[x+1,y]=0
  } else {
    X[x+1,y]=1
  }
  
  if(X[x+1,y+1]==1) {
    X[x+1,y+1]=0
  } else {
    X[x+1,y+1]=1
  }
  return(c(X))
}

equality_next_other <- function(X) # Function generating next state
                                   # Sets white on the sides if white and black in the corners
                                   # Other way round if black
{
  X <- matrix(X, sqrt(length(X)))
  
  x <- sample(2:(nrow(X)-1), 1)
  y <- sample(2:(ncol(X)-1), 1)
  
  val1 <- 0
  val2 <- 1
  if(X[x,y]==1) {
    val1 <- 1
    val2 <- 0
  }
  X[x-1,y] <- val1
  X[x,y-1] <- val1
  X[x,y+1] <- val1
  X[x+1,y] <- val1
  
  X[x-1,y-1] <- val2
  X[x-1,y+1] <- val2
  X[x+1,y-1] <- val2
  X[x+1,y+1] <- val2
  
  return(c(X))
  
}