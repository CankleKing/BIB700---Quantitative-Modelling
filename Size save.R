################################################################################
############################## Create the landscape ############################
################################################################################

### Set the seed for reproducibility, each seed represents 1 rep (N = 10)
# Set the seeds to run the model
seeds <- c(000, 111, 222, 333, 444, 555, 666, 777, 888, 999)

# Create an empty data frame to store the results
big_results <- data.frame()

# Outer loop over the seeds
for (seed in seeds) {
  
  # Inner loop over the runs
  for (run in 1:10) {
    
    # Set the seed
    set.seed(seeds)

### Create the landscape with random coordinates 
patch_name <- c("Patch.1", "Patch.2", "Patch.4", "Patch.3", "Patch.5", 
                "Patch.6", "Patch.7", "Patch.8", "Patch.9", "Patch.10")
patch_size <- as.integer(runif(10, min = 1, max = 100))
patch_coord.x <- sample(0:10, 10, replace=TRUE)
patch_coord.y <- sample(0:10, 10, replace=TRUE)
n.patch <- length(patch_name)
DatPatch <- data.frame(name=patch_name, size=patch_size, x=patch_coord.x, 
                       y=patch_coord.y)

## distance matrix
dist_patch <- as.matrix(dist(DatPatch[, c("x", "y")]))
rownames(dist_patch) <- DatPatch$name
colnames(dist_patch) <- DatPatch$name

## Plot the graph 
plot(DatPatch$x, DatPatch$y, pch=19, cex=DatPatch$size/10, col="dark green",
     xlab="X distance (Km)", ylab="Y distance (Km)")

# Print the patch names on the plot
text(DatPatch$x, DatPatch$y, labels=DatPatch$name, col="white",
     cex=DatPatch$size/100)

# add lines between patches less than 3 units away
for(i in 1:n.patch) {
  for(j in 1:n.patch) {
    if(dist_patch[i,j] < 3) {
      segments(DatPatch[i,"x"], DatPatch[i,"y"], DatPatch[j,"x"], 
               DatPatch[j,"y"], col="black", lwd = 1)
    }
  }
}
################################################################################
############################# DEFINING THE MODEL ###############################
################################################################################
############################### Local Dynamics #################################

### Create a list of data frames, with each data frame representing a patch
patch_data <- vector(mode = "list", length = n.patch) # initialize the list
for (i in 1:n.patch) {
  ## Set patch parameters
  patch_data[[i]] <- data.frame(
    N = 500, #Initial starting population size for host
    P = 50, # Initial starting population size for parasitoid
    Q = 5, # Initial starting population size for hyperparasitoid
    LAMh = (DatPatch$size[i] / 14.4), # Spatial dependent Growth rate of host 
    k1 = 0.25, # Degree of clumping for parasitoid attacks
    k2 = 1, # Degree of clumping for hyperparasitoid attacks
    kh = (0.5 * DatPatch$size[[i]]), # Host carrying capacity 
    Ap = 0.68, # Attack rate of parasitoid 
    Aq = 0.4, # Attack rate of hyperparasitoid 
    cp = 1, # Conversion efficiency of parasitoid 
    cq = 0.23, # Conversion efficiency of hyperparasitoid 
    noise = rnorm(1, mean = 1, sd = 0.1), # Stochastic noise term 
    size = DatPatch$size[[i]]) # Patch size 
}

### Define a function to update N, P, and Q
update_N_P_Q <- function(N, P, Q, dist_patch, patch_data, current_patch_) {
  LAMh <- patch_data$LAMh
  Ap <- patch_data$Ap
  Aq <- patch_data$Aq
  cp <- patch_data$cp
  cq <- patch_data$cq
  k1 <- patch_data$k1
  k2 <- patch_data$k2
  kh <- patch_data$kh
  Mu <- patch_data$Mu
  noise <- patch_data$noise
  size <- patch_data$size
  ## Calculate the new populations with local interactions 
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
  
  ## Ensure populations cannot go negative
  if (N_local < 0) N_local <- 0
  if (P_local < 0) P_local <- 0
  if (Q_local < 0) Q_local <- 0
  
  ## Return new populations as a dataframe to go to next timestep
  return(data.frame(N = N_local, P = P_local, Q = Q_local))
}

################################################################################
############################# Run the simulation ###############################
################################################################################

### Initialize a data frame to store the results
results <- data.frame(
  timestep = integer(),
  patch_name = character(),
  patch_size = integer(),
  N = integer(),
  P = integer(),
  Q = integer()
)

### Run the model for 100 timesteps
for (t in 1:100) {
  ## Loop through each patch
  for (i in 1:n.patch) {
    ## Get patch data
    patch_data_i <- patch_data[[i]]
    
    ## Update N, P, and Q
    updated_populations <- update_N_P_Q(patch_data_i$N, patch_data_i$P, 
                                        patch_data_i$Q, dist_patch, 
                                        patch_data_i)
    
    ## Update patch data
    patch_data_i$N <- updated_populations$N
    patch_data_i$P <- updated_populations$P
    patch_data_i$Q <- updated_populations$Q
    patch_data[[i]] <- patch_data_i
    
    ## Add results to data frame
    results <- rbind(results, data.frame(
      timestep = t,
      patch_name = patch_name[i],
      patch_size = patch_size[i],
      N = round(updated_populations$N),
      P = round(updated_populations$P),
      Q = round(updated_populations$Q)
    ))
  }
}

### Print the results to check data for errors 
print(results)

### Add the results table for each seed and runs to to big results table
big_results <- bind_rows(big_results, results)

  }
}


################################################################################
########################### ANALYSIS & VISUALIZATION ###########################
################################################################################

### Plotting FIG 4 - Patch size against avg. pop. size
# Load required packages
library(ggplot2)

# Calculate the average populations by patch size
avg_populations <- aggregate(cbind(N, P, Q) ~ patch_size, data = big_results, 
                             FUN = mean)

# Plot the graph
ggplot(avg_populations, aes(x = patch_size)) +
  geom_line(aes(y = N, linetype = "Host"), size = 1) +
  geom_line(aes(y = P, linetype = "Parasitoid"), size = 1) +
  geom_line(aes(y = Q, linetype = "Hyperparasitoid"), size = 1) +
  labs(x = "Patch size", y = "Population size", linetype = "Population") +
  scale_linetype_manual(values = c("Host" = "solid", "Parasitoid" = "longdash", 
                                   "Hyperparasitoid" = "dotted")) +
                        theme_classic()


### Analysis of the effect of patch size on each species population 
# Fit a GLM with Poisson distribution and log link function to populations
glm_N <- glm(N ~ patch_size, data = big_results, family = poisson(link = "log"))
glm_P <- glm(P ~ patch_size, data = big_results, family = poisson(link = "log"))
glm_Q <- glm(Q ~ patch_size, data = big_results, family = poisson(link = "log"))

# Summarize the model
summary(glm_N)
summary(glm_P)
summary(glm_Q)
