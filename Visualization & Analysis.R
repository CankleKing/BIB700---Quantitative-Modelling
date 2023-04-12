################################################################################
########################### ANALYSIS & VISUALIZATION ###########################
################################################################################
### Install required packages
install.packages("dunn.test")
install.packages("ggplot2")
install.packages("dplyr")
install.packages("dunn.test")
### Load required packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(dunn.test)


############################## Visualizing trends ##############################

### Plotting FIG 4 - Abundance over time for patches grouped by size 
## Define the ranges of patch sizes
size_ranges <- c(1, 21, 41, 61, 81, 101)

## Subset big_results to include only the first 20 time steps
big_results_subset <- big_results[big_results$timestep <= 20,]

## Bin the patches into the size ranges
big_results_subset$size_range <- cut(big_results_subset$patch_size, 
                                     breaks = size_ranges,right = FALSE)

## Calculate mean and standard deviation of population sizes for each size range
host_size_means_N <- aggregate(N ~ timestep + size_range, 
                               data = big_results_subset, FUN = mean) 
host_size_sds_N <- aggregate(N ~ timestep + size_range, 
                             data = big_results_subset, FUN = sd)
host_size_means_P <- aggregate(P ~ timestep + size_range, 
                               data = big_results_subset, FUN = mean) 
host_size_sds_P <- aggregate(P ~ timestep + size_range, 
                             data = big_results_subset, FUN = sd)
host_size_means_Q <- aggregate(Q ~ timestep + size_range, 
                               data = big_results_subset, FUN = mean) 
host_size_sds_Q <- aggregate(Q ~ timestep + size_range, 
                             data = big_results_subset, FUN = sd)

## Set theme for all graphs
my_theme <- theme_bw() + 
  theme(panel.grid.major = element_line(colour = "gray90", linetype = "dashed"),
        panel.grid.minor = element_line(colour = "gray95", linetype = "dashed"))

### Plot the mean host population size over time for each size range (FIG 4a)
ggplot(host_size_means_N, aes(x = timestep, y = N, fill = size_range)) + 
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(y = N, ymin = N - host_size_sds_N$N, 
                    ymax = N + host_size_sds_N$N), 
                width = 0.2, position = position_dodge(0.9)) +
  scale_fill_brewer(palette = "Set3") +
  labs(x = "Time", y = "Mean U. cardui Population Size") + 
  my_theme
### Plot mean parasitoid population size over time for each size range (FIG 4b)
ggplot(host_size_means_P, aes(x = timestep, y = P, fill = size_range)) + 
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(y = P, ymin = P - host_size_sds_P$P, 
                    ymax = P + host_size_sds_P$P), 
                width = 0.2, position = position_dodge(0.9)) +
  scale_fill_brewer(palette = "Set3") +
  labs(x = "Time", 
       y = "Mean E. serratulae Population Size") +
  my_theme
### Plot mean hyperparasitoid pop. size over time for each size range (FIG 4c)
ggplot(host_size_means_Q, aes(x = timestep, y = Q, fill = size_range)) + 
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(y = Q, ymin = Q - host_size_sds_Q$Q, 
                    ymax = Q + host_size_sds_Q$Q), 
                width = 0.2, position = position_dodge(0.9)) +
  scale_fill_brewer(palette = "Set3") +
  labs(x = "Time", 
       y = "Mean E. robusta Population Size") +
  my_theme



################################## Summary table ###############################
## Bin the patches into the size ranges
big_results$size_range <- cut(big_results$patch_size, breaks = size_ranges, 
                              right = FALSE)

### Group the dataset by size range and compute summary statistics
summary_table <- big_results %>% 
  group_by(size_range) %>% 
  summarize(mean_N = round(mean(N), 2),
            sd_N = round(sd(N), 2),
            min_N = min(N),
            max_N = max(N),
            mean_P = round(mean(P), 2),
            sd_P = round(sd(P), 2),
            min_P = min(P),
            max_P = max(P),
            mean_Q = round(mean(Q), 2),
            sd_Q = round(sd(Q), 2),
            min_Q = min(Q),
            max_Q = max(Q))
## Print the summary table - TABLE 2 
summary_table



################################# Kruskal-wallis ###############################
## Calculate the average populations by patch size
avg_population <- aggregate(cbind(N, P, Q) ~ patch_size, 
                            data = big_results,
                            FUN = mean)

## Bin the patches into the size ranges
avg_population$size_range <- cut(avg_population$patch_size, 
                                 breaks = size_ranges, 
                                 right = FALSE)

### Check the assumptions for One-way ANOVA
## 1. Normality assumption
shapiro.test(avg_population$N) # output P-value > 0.05 so normally distributed 
shapiro.test(avg_population$P) # output P-value < 0.05 so non-parametric needed
shapiro.test(avg_population$Q) # output P-value < 0.05 so non-parametric needed

## 2. Homogeneity of variance assumption
bartlett.test(N ~ patch_size, data = big_results_t20up) # Assumption violated

# Therefore we do a non-parametric kruskall-wallis 
## Perform Kruskal-Wallis test
K_N_results <- kruskal.test(N ~ size_range, data = avg_population)
K_P_results <- kruskal.test(P ~ size_range, data = avg_population)
K_Q_results <- kruskal.test(Q ~ size_range, data = avg_population)

## Perform Dunn post-hoc test
dunn_resultsN <- dunn.test(avg_population$N, avg_population$size_range, 
                            method = "bonferroni")
dunn_resultsP <- dunn.test(avg_population$P, avg_population$size_range,
                            method = "bonferroni")
dunn_resultsQ <- dunn.test(avg_population$Q, avg_population$size_range,
                            method = "bonferroni")



########################## Simple Linear regression ############################
### Plotting FIG 5 - Patch size against avg. pop. size

## Subset data for average populations w/out initial population fluctuations
big_results_t20up <- big_results[big_results$timestep <= 20 &
                                   big_results$patch_size != 61,] # - outlier

## Calculate the average populations by patch size
avg_population <- aggregate(cbind(N, P, Q) ~ patch_size, 
                            data = big_results_t20up,
                            FUN = mean)

## Fitting linear regression models to the data
model_N <- lm(N ~ patch_size, data = avg_population)
summary(model_N) # View the model summary

model_P <- lm(P ~ patch_size, data = avg_population)
summary(model_P) # View the model summary

model_Q <- lm(Q ~ patch_size, data = avg_population)
summary(model_Q) # View the model summary

### Checking the assumptions
## Linearity
## Visual check of linearity for each pop from scatter & residual vs fitted plot
# Hosts
ggplot(avg_population, aes(x = patch_size, y = N)) +
  geom_point() + # Add the points
  geom_smooth(method = "lm") + # Add the linear regression line
  labs(x = "Patch size", y = "Average N") # Add axis labels
plot(model_N,1)
# Parasitoids
ggplot(avg_population, aes(x = patch_size, y = P)) +
  geom_point() + 
  geom_smooth(method = "lm") + 
  labs(x = "Patch size", y = "Average P") 
plot(model_P,1)
# Hyperparasitoids
ggplot(avg_population, aes(x = patch_size, y = P)) +
  geom_point() + 
  geom_smooth(method = "lm") + 
  labs(x = "Patch size", y = "Average Q") 
plot(model_Q,1)

## Homogeneity of variance 
plot(model_N, 3)
plot(model_P, 3)
plot(model_Q, 3) # Normal distribution for all so assumption is met

## Normality of residuals 
plot(model_N, 2)
plot(model_P, 2)
plot(model_Q, 2) # Normality assumption is met 

## Check for outliers and high leverage points
plot(model_N, 5)
plot(model_P, 5)
plot(model_Q, 5) # The data point for patch size = 62 has a high leverage 
# Therefore it can be removed from the dataset (This has already been done)

## Extracting R-squared values for each model
rsq_N <- summary(model_N)$r.squared
rsq_P <- summary(model_P)$r.squared
rsq_Q <- summary(model_Q)$r.squared

### Plotting all three scatter plots and linear regressions on the same graph
ggplot(avg_population, aes(x = patch_size)) +
  geom_point(aes(y = N, color = "N"), size = 3, shape = "+") +
  geom_smooth(aes(y = N), method = "lm", se = FALSE, color = "black",
              lty = "solid") +
  geom_point(aes(y = P, color = "P"), size = 3, shape = "+") +
  geom_smooth(aes(y = P), method = "lm", se = FALSE, color = "blue", 
              lty = "solid") +
  geom_point(aes(y = Q, color = "Q"), size = 3, shape = "+") +
  geom_smooth(aes(y = Q), method = "lm", se = FALSE, color = "red", 
              lty = "solid") +
  labs(x = "Patch size", y = "Population size", color = "Population") +
  scale_color_manual(values = c("N" = "black", "P" = "blue", "Q" = "red")) +
  my_theme +
  geom_text(aes(x = 20, y = 35, label = (paste0("R^2 = ", 
                                                round(rsq_N, 3),
                                                ";   Pr(>|t|) = < 2.2e-16"))),
            color = "black", size = 3, vjust = -1) +
  geom_text(aes(x = 20, y = 30, label = (paste0( "R^2 = ", 
                                                round(rsq_P, 3),
                                                ";   Pr(>|t|) = < 2.2e-16"))),
            color = "blue", size = 3, vjust = -1) +
  geom_text(aes(x = 20, y = 25, label = (paste0("R^2 = ", 
                                                round(rsq_Q, 3),
                                                ";   Pr(>|t|) = < 2.2e-16"))),
            color = "red", size = 3, vjust = -1)



######################### Generalised Linear Model #############################

### Fit a GLM to populations with all possible interactions 
glm_N <- glm(N ~ patch_size + P + Q + P*Q + patch_size*P + patch_size*Q, 
             data = avg_population, family = poisson(link = "log"))

glm_P <- glm(P ~ patch_size + N + Q + N*Q + patch_size*N + patch_size*Q, 
             data = avg_population, family = poisson(link = "log"))

glm_Q <- glm(Q ~ patch_size + P + N + N*P + patch_size*P + patch_size*N, 
             data = avg_population, family = poisson(link = "log"))

## Summarize the model
summary(glm_N)
summary(glm_P)
summary(glm_Q)

## Check data for model for over dispersion (supplement in each model)
residual_deviance <- summary(glm_Q)$deviance 
residual_df <- summary(glm_Q)$df.residual 
dispersion <- residual_deviance / residual_df
dispersion
# none of the models were overdispersed 

