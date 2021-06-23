# This script contains the simulation methods
# First the different methods of simulating data built, then they are used for the actual simulation

rNegBin <- function(otu=NULL, n)
{
  ##############################################################################
  # otu:
  #   Input column from an existing otu table
  # n:
  #   Number of draws from the negative binomial distribution
  ##############################################################################
  # Output: 
  #   n counts taken from a single variate negative binomial distribution.
  ##############################################################################
  
  otuMean <- mean(otu)
  otuVar <- var(otu)
  
  otuSize <- otuMean**2 / (otuVar - otuMean)
  if(otuSize <= 0){stop("Size parameter smaller than or equal to zero.")}
  
  rnbinom(n, mu = otuMean, size = otuSize)
}

simData <- function(simMethods, n, N)
{
  ##############################################################################
  # simMethods:
  #   A vector containing strings depicting what simulation methods to use.
  # n:
  #   Number of draws per otu for each simulation.
  # N:
  #   Number of simulations.
  ##############################################################################
  # Output: 
  #   A list containing the simulated data.
  ##############################################################################
  
}