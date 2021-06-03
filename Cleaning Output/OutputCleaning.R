# Cleaning simulation output
library(tidyverse)
library(ggpubr)

##### Covariance structure simulation output #####

# Setting up a custom color scale for the line colors
cc <- scales::seq_gradient_pal("white","black", "Lab")(seq(.1,1,length.out = 8))

# Function to load and tidy covariance structure simulation results
tidy_dimension_rda <- function(f){
  df <- readRDS(f)
  df_tidy <- data.frame(t(df[,-1]))
  df_tidy <- cbind(1:10, df_tidy)
  colnames(df_tidy) <- c('dim', df[, 1])
  df_tidy <- pivot_longer(df_tidy, !dim, names_to = 'rho', values_to = 'val')
}

# Load simulation results (in this case, the compound symmetric simulation results)
cs_means <- tidy_dimension_rda("cs_means.Rda")
cs_vars <- tidy_dimension_rda("cs_vars.Rda")
cs_rss <- tidy_dimension_rda("cs_rss.Rda")
cs_test_rss<- tidy_dimension_rda("cs_test_rss.Rda")

# Combine loaded results into a single data frame for plotting
cs_tidy <- data.frame(dim = cs_means$dim, rho = cs_means$rho, 
                             edf = cs_means$val, sd = (cs_vars$val)^(1/2), 
                             rss = cs_rss$val, test_rss = cs_test_rss$val)
# Plot EDF with 1-sd bars
ggplot(cs_tidy, aes(x = dim, color = rho, linetype = rho)) + 
  geom_line(aes(y = edf)) +
  geom_errorbar(aes(ymin = edf -sd, ymax = edf + sd), position=position_dodge(0.1), width=.4) +
  labs(title = 'Compound Symmetric Case', 
       x = "Envelope Dimension", y = "EDF", color = "\u03C1", linetype = "\u03C1") +
  scale_x_continuous(breaks = seq(1,10)) +
  theme_minimal(base_size = 12) + 
  theme(legend.title=element_text(size=14, face = "bold"), 
        legend.text=element_text(size=14)) + 
  scale_color_manual(values = cc) +
  scale_linetype_manual(values = rep(c("solid", "dashed", "dotted"), 3))

# Plot EDF with RSS and in-sample test RSS (for rho = 0.8)
cs_tidy_sample <- filter(cs_tidy, rho == 0.8)

ggplot(cs_tidy_sample, aes(x = dim)) + 
  geom_line(aes(y = edf, linetype = "EDF")) +
  geom_line(aes(y = rss, linetype = "RSS")) +
  geom_line(aes(y = test_rss, linetype = "Test RSS")) +
  labs(title = 'Envelope, RSS, and Test RSS in Compound Symmetric Case', 
       x = "Envelope Dimension", y = "") +
  scale_x_continuous(breaks = seq(1,10)) +
  theme_minimal(base_size = 12) + 
  theme(legend.title=element_blank(), 
        legend.text=element_text(size=14))

##### Eigenvalue simulation output #####
# Setting up a custom color scale for the line colors
cc2 <- scales::seq_gradient_pal("white","black", "Lab")(seq(.3,1,length.out = 6))

# Function to load and tidy eigenvalue simulation results
tidy_eigen_rda <- function(f){
  df <- readRDS(f)
  df_tidy <- data.frame(t(df[,-1]))
  df_tidy <- cbind(1:10, df_tidy)
  colnames(df_tidy) <- c('dim',df[,1])
  df_tidy <- pivot_longer(df_tidy, !dim, names_to = 'first_eigen', values_to = 'val')
}

# Load simulation results
eigen_means <- tidy_eigen_rda("eigen_means.Rda")
eigen_vars <- tidy_eigen_rda("eigen_vars.Rda")
eigen_rss <- tidy_eigen_rda("eigen_rss.Rda")
eigen_test_rss <- tidy_eigen_rda("eigen_test_rss.Rda")

# Combine loaded results into a single data frame for plotting
eigen_tidy <- data.frame(dim = eigen_means$dim, 
                         first_eigen = factor(eigen_means$first_eigen, levels = c('1', '5', '10', '25', '100') ), 
                         edf = eigen_means$val, sd = (eigen_vars$val)^(1/2),
                         rss = eigen_rss$val, test_rss = eigen_test_rss$val)

# Plot EDF with 1-sd bars
ggplot(eigen_tidy, aes(x = dim, y = edf, color = first_eigen, linetype = first_eigen)) + 
  geom_line() +
  geom_text( aes(label = ifelse( round(edf,1) < 10, round(edf,1), '') ), hjust=-1, vjust=0, show.legend = F) +
  geom_errorbar(aes(ymin = edf-sd, ymax = edf + sd, color = first_eigen), position=position_dodge(0.1), width=.4) +
  labs(title = "Reduction Case 1", x = "Envelope Dimension", y = "EDF", 
       color = expression( paste(lambda[1]) ), 
       linetype =  expression( paste(lambda[1]) ) 
       ) +
  scale_x_continuous(breaks = seq(1,10)) +
  theme_minimal(base_size = 12) +
  theme(legend.title=element_text(size=14, face = "bold"), 
        legend.text=element_text(size=14),
        legend.title.align = 0.5) +
  scale_color_manual(values = cc2) +
  scale_linetype_manual(values = rep(c("solid", "dashed", "dotted"), 2))

# Plot EDF with RSS and in-sample test RSS (for lambda1 = 25)
eigen_tidy_sample <- filter(eigen_tidy, first_eigen == 25)

ggplot(eigen_tidy_sample, aes(x = dim)) + 
  geom_line(aes(y = edf, linetype = "EDF")) +
  geom_line(aes(y = rss, linetype = "RSS")) +
  geom_line(aes(y = test_rss, linetype = "Test RSS")) +
  labs(title = 'Envelope, RSS, and Test RSS in Reduction Case 1', 
       x = "Envelope Dimension", y = "") +
  scale_x_continuous(breaks = seq(1,10)) +
  theme_minimal(base_size = 12) + 
  theme(legend.title=element_blank(), 
        legend.text=element_text(size=14))

##### Minneapolis simulation output #####

# Load data from Minneapolis predictors only simulation
mpls_means <- readRDS("mpls_means.Rda")
mpls_vars <- readRDS("mpls_vars.Rda")
mpls_tidy <- data.frame(dim= 1:5, edf = mpls_means, sd = (mpls_vars)^(1/2)  )

# Plot EDF with 1-sd bars
ggplot(mpls_tidy, aes(dim, edf)) + geom_line( size = 0.7) +
  geom_errorbar(aes(ymin = edf -sd, ymax = edf + sd), position=position_dodge(0.1), width=.1, linetype = 2) +
  labs(#title = 'Envelope EDF with Real Covariance Matrix', 
    x = "Envelope Dimension", y = "EDF") +
  scale_x_continuous(breaks = seq(1,5)) +
  scale_y_continuous(breaks = seq(1,12)) +
  coord_cartesian(ylim=c(0,12)) +
  theme_minimal(base_size = 14)

# Load data from Minneapolis parametric bootstrap simulation
mpls_boot_means <- readRDS("mpls_bootstrap_means.Rda")
mpls_boot_tidy <- data.frame(dim= 1:5, edf = mpls_boot_means)

# Plot EDF with 1-sd bars
ggplot(mpls_boot_tidy, aes(dim, edf)) + geom_line( size = 0.7) +
  labs(#title = 'Envelope EDF with Real Covariance Matrix', 
    x = "Envelope Dimension", y = "EDF") +
  scale_x_continuous(breaks = seq(1,5)) +
  scale_y_continuous(breaks = seq(1,12)) +
  coord_cartesian(ylim=c(0,12)) +
  theme_minimal(base_size = 14)
