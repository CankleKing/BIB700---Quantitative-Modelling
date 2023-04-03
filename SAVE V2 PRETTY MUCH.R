############################## Create the landscape ############################

# Create the landscape with random coordinates 
patch_name <- c("Patch.7", "Patch.6", "Patch.5", "Patch.3", "Patch.4", 
                "Patch.2", "Patch.1")
patch_size <- sample(30:100, 7, replace=TRUE)
patch_coord.x <- sample(0:10, 7, replace=TRUE)
patch_coord.y <- sample(0:10, 7, replace=TRUE)
n.patch <- length(patch_name)
DatPatch <- data.frame(name=patch_name, size=patch_size, x=patch_coord.x, 
                       y=patch_coord.y)

# distance matrix
dist_patch <- as.matrix(dist(DatPatch[, c("x", "y")]))
rownames(dist_patch) <- DatPatch$name
colnames(dist_patch) <- DatPatch$name


plot(DatPatch$x, DatPatch$y, pch=19, cex=DatPatch$size/10, col="dark green")

# Print the patch names on the plot
text(DatPatch$x, DatPatch$y, labels=DatPatch$name, col="white",
     cex=DatPatch$size/100)

# add lines between patches less than 2.4 units away
for(i in 1:n.patch) {
  for(j in 1:n.patch) {
    if(dist_patch[i,j] < 3) {
      segments(DatPatch[i,"x"], DatPatch[i,"y"], DatPatch[j,"x"], 
               DatPatch[j,"y"], col="black", lwd = 1)
    }
  }
}

############################# Local Dynamics stage #############################

### Create a list of data frames, with each data frame representing a patch
patch_data <- vector(mode = "list", length = n.patch) # initialize the list
for (i in 1:n.patch) {
  ## Set patch parameters
  patch_data[[i]] <- data.frame(
    N = 500, #Initial starting population size for host
    P = 50, # Initial starting population size for parasitoid
    Q = 5, # Initial starting population size for hyperparasitoid
    LAMh = (DatPatch$size[i] / 5.6), # Spatial dependent Growth rate of host 
    k1 = 0.25, # Degree of clumping for parasitoid attacks
    k2 = 1, # Degree of clumping for hyperparasitoid attacks
    kh = (5.6 * DatPatch$size[[i]]), # Host carrying capacity 
    Ap = 0.6, # Attack rate of parasitoid 
    Aq = 0.4, # Attack rate of hyperparasitoid 
    cp = 1, # Conversion efficiency of parasitoid 
    cq = 0.23, # Conversion efficiency of hyperparasitoid 
    noise = rnorm(1, mean = 1, sd = 0.2), # Stochastic noise term 
    Mu = 0.0001,
    ## Include patch size and connectivity for later 
    conn = length(DatPatch$name[which(dist_patch[[i]] < 3)]), # No. of neighbors
    size = DatPatch$size[[i]])
}

### Define a function to update N, P, and Q
update_N_P_Q <- function(N, P, Q, dist_patch, patch_data) {
  LAMh <- patch_data$LAMh
  Ap <- patch_data$Ap
  Aq <- patch_data$Aq
  cp <- patch_data$cp
  cq <- patch_data$cq
  k1 <- patch_data$k1
  k2 <- patch_data$k2
  kh <- patch_data$kh
  Mu <- patch_data$Mu
  noise <-patch_data$noise
  ## Calculate the new popualtions with local interactions 
  Host_dens_dep <- LAMh * exp((-N * exp(-Ap*P)) / kh)
  N_local <- Host_dens_dep * N * ((1 + ((Ap * P) / k1)) ^ - k1) * 
    ((1 + ((Aq * Q) / k2)) ^ - k2) * noise
  P_local <- cp * N * ((1 + ((Aq * Q) / k2)) ^ - k2) * 
    ((1 - ((1 + ((Ap * P) / k1)) ^ - k1))) * noise
  Q_local <- cq * N * ((1 - (1 + ((Aq * Q) / k2)) ^ - k2)) * noise
  
  ## Set a cap on the population size 
  N_local <- min(N_local, kh * noise)
  P_local <- min(P_local, kh * noise)
  Q_local <- min(Q_local, kh * noise)
  
  return(data.frame(N = N_local, P = P_local, Q = Q_local))
  
}

################################# Dispersal Stage ##############################

### Create a loop of for and if statements to model dispersal between patches 
for (i in 1:n.patch) {
  ## get the names of connected patches
  connected_patches <- DatPatch$name[which(dist_patch[i,] < 3)]
  
### if there are connected patches, then disperse some of the population
  if (length(connected_patches) > 0) {
    ## calculate the amount to disperse based on Mu
    N_disperse <- patch_data[[i]]$N * patch_data[[i]]$Mu
    P_disperse <- patch_data[[i]]$P * patch_data[[i]]$Mu
    Q_disperse <- patch_data[[i]]$Q * patch_data[[i]]$Mu
    
    ## randomly allocate some of the population to each connected patch
    N_dispersed <- sample(c(N_disperse, rep(0, length(connected_patches)-1)))
    P_dispersed <- sample(c(P_disperse, rep(0, length(connected_patches)-1)))
    Q_dispersed <- sample(c(Q_disperse, rep(0, length(connected_patches)-1)))
    total_dispersed <- N_dispersed + P_dispersed + Q_dispersed
    proportion_dispersed <- total_dispersed / sum(total_dispersed)
    
    ## update the population of each connected patch
    for (j in 1:length(connected_patches)) {
      patch_data[[which(DatPatch$name == connected_patches[j])]]$N <- 
        patch_data[[which(DatPatch$name == connected_patches[j])]]$N 
      + N_dispersed[j] * proportion_dispersed[j]
      patch_data[[which(DatPatch$name == connected_patches[j])]]$P <- 
        patch_data[[which(DatPatch$name == connected_patches[j])]]$P 
      + P_dispersed[j] * proportion_dispersed[j]
      patch_data[[which(DatPatch$name == connected_patches[j])]]$Q <- 
        patch_data[[which(DatPatch$name == connected_patches[j])]]$Q 
      + Q_dispersed[j] * proportion_dispersed[j]
    }
    
    ## subtract the amount dispersed from the current patch
    patch_data[[i]]$N <- patch_data[[i]]$N - N_disperse
    patch_data[[i]]$P <- patch_data[[i]]$P - P_disperse
    patch_data[[i]]$Q <- patch_data[[i]]$Q - Q_disperse
  }
}

# Run the model for each patch
num_time_steps <- 100
for (i in 1:num_time_steps) {
  for (j in 1:n.patch) {
    # Update the populations of the current patch
    N_P_Q <- update_N_P_Q(patch_data[[j]]$N, patch_data[[j]]$P, 
                          patch_data[[j]]$Q, dist_patch, patch_data[[j]])
    patch_data[[j]][c("N", "P", "Q")] <- N_P_Q
    print(paste("Patch", DatPatch$name[j], "N =", round(N_P_Q$N), "; P =", 
                round(N_P_Q$P), "; Q =", round(N_P_Q$Q)))
  }
}


## Create an empty table for the parameter variables 
Results <- data.frame(Time_Step = t,
                      Patch_Name = patch_data$name[[i]],
                      Patch_size = patch_data[[i]],
                      Connectivity = patch_data[[i]]$conn,
                      Dispersal = patch_data[[i]]$Mu,
                      N_Abundance = patch_data[[i]]$N,
                      P_Abundance = patch_data[[i]]$P,
                      Q_Abundance = patch_data[[i]]$Q,
                      stringsAsFactors = FALSE)

