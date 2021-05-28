# Code to simulate data for the analysis -------------------------------
library(tidyverse)
library(readr)
library(mice)

# The aim of this script is to simulate the synthetic dataset that will
# be used as basis to reproduce the workflow of the analysis implemented
# in the study.
set.seed(745698)

# 1) Generate baseline variables ---------------------------------------
# 30 baseline variables will be simulated from a multivariate normal
# distribution. 15 of them will be dichotomized to obtain binary
# categorical variables.
n <- 500
p <- 14

mus <- rep(0, p)
sds <- rep(1, p)
corr_mat <- matrix(rep(0, p * p), nrow = p, ncol = p)
diag(corr_mat) <- 1
corrs <- runif(n = sum(lower.tri(corr_mat)), min = 0.1, max = 0.3)

corr_mat[upper.tri(corr_mat)] <- corrs
corr_mat[lower.tri(corr_mat)] <- t(corr_mat)[lower.tri(t(corr_mat))]
corr_matrix <- as.matrix(Matrix::nearPD(corr_mat)$mat)
cov_mat <- sweep(sweep(corr_matrix, 1, sds, "*"), 2, sds, "*")

# Generate X matrix
x_mat <- MASS::mvrnorm(n = n, mu = mus, Sigma = cov_mat)

# Make half of X binary
x_mat[, ((p/2) + 1):p] <- apply(
   x_mat[, ((p/2) + 1):p], 2, function(x) ifelse(x < 0, 0, 1)
)

# Add the intercept and return the design matrix
design_matrix <- cbind(rep(1, n), x_mat)

# 2) Simulate the treatment assignment ---------------------------------
# Simulate the treatment assignment according to a logistic regression
# model with main effects only. Coefficients will be simulated from
# uniform distribution
betas <- c(-1.8, runif(n = p, min = log(0.8), max = log(1.75)))
lp <- design_matrix %*% betas
prob <- exp(lp)/(1 + exp(lp))
treat <- stats::rbinom(n = n, size = 1, prob = prob)

# Make final tibble with complete observed data ------------------------
d_complete <- design_matrix %>%
   tibble::as_tibble() %>%
   dplyr::mutate(treat = treat) %>%
   dplyr::mutate(across(c(V9:V15), ~ factor(.))) %>%
   dplyr::select(-V1)

# 3) Simulate missing data with ampute function ------------------------
# Initialize the function ----------------------------------------------
init <- mice::ampute(d_complete %>% dplyr::select(-treat))

# Define the mechanism and the proportion of missingness
mechanism <- "MAR"
prop <- 0.15

# Simulate 7 patterns of missingness -----------------------------------
n_patt <- 7
patterns <- matrix(
   n_patt %>%
      purrr::rerun(stats::rbinom(n = p, size = 1, prob = 0.5)) %>%
      unlist(),
   ncol = p,
   nrow = n_patt
)

# Define the frequencies of each pattern ------------------------------
freqs <- rep(1/n_patt, n_patt)

# Simulate the weights for each variable in the patterns
sim_wts <- stats::runif(
   length(patterns != 0),
   min = log(1.01),
   max = log(1.8)
)

wts <- ifelse(patterns == 0, 0, sim_wts)

# Generate missing data -----------------------------------------------
d_incomplete <- mice::ampute(
   data = d_complete %>% dplyr::select(-treat),
   patterns = patterns,
   freq = freqs,
   weights = wts,
   prop = prop,
   mech = mechanism,
   cont = TRUE,
   type = "RIGHT"
)

# Create the final dataset --------------------------------------------
d <- d_incomplete$amp %>%
   dplyr::mutate(treat = d_complete$treat) %>%
   dplyr::mutate(
      dplyr::across(
         c(V9:V15),
         ~ factor(dplyr::if_else(. == 1, "0", "1"), levels = c("0", "1"))
      )
   )

# 4) Save the dataset as rds ------------------------------------------
readr::write_rds(x = d, file = here::here("data", "synthetic_data.rds"))

