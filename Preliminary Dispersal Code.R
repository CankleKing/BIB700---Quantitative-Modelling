################################################################################
##################### Preliminary dispersal methodology ########################
################################################################################

### Create the random landscape 
patch_name <- c("Patch.7", "Patch.6", "Patch.5", "Patch.3", "Patch.4", "Patch.2"
                , "Patch.1") # Name patches
patch_size <- sample(30:100, 7, replace=TRUE) # Set range of patch sizes 
patch_coord.x <- sample(0:10, 7, replace=TRUE) # Set range of x coordinates
patch_coord.y <- sample(0:10, 7, replace=TRUE) # Set range of y coordinates
n.patch <- length(patch_name)
DatPatch <- data.frame(name=patch_name, size=patch_size, x=patch_coord.x, 
                       y=patch_coord.y)

### Create a distance matrix to work out connected patches
dist_patch <- as.matrix(dist(DatPatch[, c("x", "y")]))
rownames(dist_patch) <- DatPatch$name
colnames(dist_patch) <- DatPatch$name

### Plot the landscape
plot(DatPatch$x, DatPatch$y, pch=19, cex=DatPatch$size/10, col="dark green")

### Print the patch names on the plot
text(DatPatch$x, DatPatch$y, labels=DatPatch$name, col="white",
     cex=DatPatch$size/100)

### Add lines between patches less than 3km away
for(i in 1:n.patch) {
  for(j in 1:n.patch) {
    if(dist_patch[i,j] < 3) {
      segments(DatPatch[i,"x"], DatPatch[i,"y"], DatPatch[j,"x"], 
               DatPatch[j,"y"], col="black", lwd = 1)
    }
  }
}

### Create a list of data frames, with each data frame representing a patch
patch_data <- vector(mode = "list", length = n.patch) # initialize the list
for (i in 1:n.patch) {
  ## set patch parameters for each patch 
  patch_data[[i]] <- data.frame(
    N = 500, #Initial starting population size for host
    P = 50, # Initial starting population size for parasitoid
    Q = 5, # Initial starting population size for hyperparasitoid
    LAMh = (DatPatch$size[i] / 14.4), # Spatialy dependent Growth rate of host 
    k1 = 0.25, # Degree of clumping for parasitoid attacks
    k2 = 1, # Degree of clumping for hyperparasitoid attacks
    kh = (0.5 * DatPatch$size[i]), # Spatially dependent host carrying capacity 
    Ap = 0.6, # Attack rate of parasitoid 
    Aq = 0.4, # Attack rate of hyperparasitoid 
    cp = 1, # Conversion efficiency of parasitoid 
    cq = 0.23, # Conversion efficiency of hyperparasitoid 
    noise = rnorm(1, mean = 1, sd = 0.2), # Stochastic noise term 
    MuN = 0.1, # Dispersal fraction of hosts  
    MuP = 0.1, # Dispersal fraction of parasitoids 
    MuQ = 0.1, # Dispersal fraction of hyperparasitoids 
    
    ## Include patch size and connectivity for later 
    conn = length(DatPatch$name[which(dist_patch[[i]] < 3)]), # No. of neighbors
    size = DatPatch$size[i]
  )
}


############## Get the average populations from connected patches ##############
# Initialize empty data frame
patch_neighbours <- data.frame(patch=character(), neighbour_patch=character(),
                               distance=numeric(), N=numeric(), P=numeric(), 
                               Q=numeric())

# Loop over each patch
for (i in 1:n.patch) {
  # Get the neighbors within 3 units or less
  neighbours <- DatPatch$name[which(dist_patch[i,] < 3)]
  neighbours <- neighbours[neighbours != DatPatch$name[i]] # Remove itself
  
  # Initialize empty data frame to store the neighbor information
  neighbour_info <- data.frame(patch=character(), connected=character(), 
                               N=numeric(), P=numeric(), Q=numeric())
  
  # Loop over the neighbors and get their populations of N, P, Q
  for (j in neighbours) {
    neighbour_data <- patch_data[[which(DatPatch$name == j)]]
    neighbour_info <- rbind(neighbour_info, data.frame(connected=j, 
                                                       N=neighbour_data$N, 
                                                       P=neighbour_data$P, 
                                                       Q=neighbour_data$Q))  
    
    # Add the current patch's neighbor information to the data frame
    patch_neighbours <- rbind(patch_neighbours, 
                              data.frame(patch=DatPatch$name[[i]], 
                                         neighbour_patch=j, N=neighbour_info$N, 
                                         P=neighbour_info$P, Q=neighbour_info$Q))
  }
  
  # Check if the patch has no neighbors
  if (length(neighbours) == 0) {
    patch_neighbours <- rbind(patch_neighbours, 
                              data.frame(patch=DatPatch$name[[i]], 
                                         neighbour_patch=NA, N=0, P=0, Q=0))
  }
}

# Print the resulting data frame to check for errors
patch_neighbours

# Sum all the connected patches
patch_neighbours_sum <- aggregate(cbind(N, P, Q) ~ patch, 
                                  data = patch_neighbours, FUN = sum)


# Create an empty data frame with = same number of rows as 'patch_neighbours_sum'
av_N_P_Q <- data.frame(matrix(0, nrow = nrow(patch_neighbours_sum), ncol = 3))

# Name the columns of the new data frame
colnames(av_N_P_Q) <- c("av_N", "av_P", "av_Q")

# Calculate the average of all connected patches for each patch
for (i in 1:nrow(patch_neighbours_sum)) {
  av_N_P_Q[i, "av_N"] <- patch_neighbours_sum[i, "N"] / length(unique(c(patch_neighbours$neighbour_patch[patch_neighbours$patch == patch_neighbours_sum[i, "patch"]], patch_neighbours$patch[patch_neighbours$neighbour_patch == patch_neighbours_sum[i, "patch"]])))
  av_N_P_Q[i, "av_P"] <- patch_neighbours_sum[i, "P"] / length(unique(c(patch_neighbours$neighbour_patch[patch_neighbours$patch == patch_neighbours_sum[i, "patch"]], patch_neighbours$patch[patch_neighbours$neighbour_patch == patch_neighbours_sum[i, "patch"]])))
  av_N_P_Q[i, "av_Q"] <- patch_neighbours_sum[i, "Q"] / length(unique(c(patch_neighbours$neighbour_patch[patch_neighbours$patch == patch_neighbours_sum[i, "patch"]], patch_neighbours$patch[patch_neighbours$neighbour_patch == patch_neighbours_sum[i, "patch"]])))
}

# Add the new columns to the 'patch_neighbors_sum' data frame
patch_neighbours_sum <- cbind(patch_neighbours_sum, av_N_P_Q)

################################################################################
# The final version of the table patch_neighbours_sum contains the values of the
# average populations from all connected patches.
# If ye' mighty adventurer feel brave enough please venture fourth and attempt 
# to defeat the mighty dragon that is: coding a dispersal phase using updated
# populations after local dynamics and then running this through each timestep
# If you manage please send me a message with your method!!!
################################################################################


### Define a function to update N, P, and Q for dispersal stage
update_N_P_Q2 <- function(N, P, Q, dist_patch, patch_data, current_patch_idx) {
  ## Loads parameters with update to avoid errors 
  MuN <- patch_data$MuN
  MuP <- patch_data$MuP
  MuQ <- patch_data$MuQ
  noise <-patch_data$noise
  conn <- patch_data$name
  size <- patch_data$size
  
  ### calculate the average population across connected patches
  
  AvN <- patch_neighbours_sum$av_N
  AvP <- patch_neighbours_sum$av_P
  AvQ <- patch_neighbours_sum$av_Q
  ## calculate the dispersed populations: Equation No. 5 - Dispersal 
  N_disp <- (((1 - MuN) * patch_data$N) + (MuN * AvN))
  P_disp <- (((1 - MuP) * patch_data$P) + (MuP * AvP))
  Q_disp <- (((1 - MuQ) * patch_data$Q) + (MuQ * AvQ))
  
  return(data.frame(N = N_disp, P = P_disp, Q = Q_disp))
  
}
