# This is the main script of the benchmark architecture
# Test version 23/06/2021

# Load packages
library(phyloseq)

# Benchmark settings:
#   Simulation settings:
minLibSize <- 100
test <- TRUE
simMethods <- c("negative binomial no correlation")
N <- 10 # number of simulations
n <- 100 # number of samples per simulation
phyloseqData <- NULL # initialization variable name of the data object, will get an actual value when loading data
simulatedData <- NULL  # initialization variable name of the simulated data object, will get an actual value when simulating data

# Load data
if (test)
  source("load_testData.R")
  
# Preprocess data
source("data_prepocessing.R")
phyloseqData <- preprocData(phyloseqData, minLibSize)

# Simulate data
source("simulate_data.R")
simulatedData <- simData(simMethods=simMethods, n=n, N=N)





