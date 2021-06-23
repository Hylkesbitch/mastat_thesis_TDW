# This script is used to load the datasets

# To develop the benchmark the HMP16S data is used
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

# The following initializes usage of Bioc devel
#BiocManager::install(version = '3.12')

#BiocManager::install("HMP16SData")
#BiocManager::install("phyloseq")
library(phyloseq)
library(HMP16SData)
V13_phyloseq <- as_phyloseq(V13())
libSizes <- colSums(V13_phyloseq@otu_table)

#-------------------------------------------------------------------------------
# PARAMETRIC SIMULATION
## NEGATIVE BINOMIAL

# SpiecEasi is a package designed to infer graphical models for microbiome 
# relative abundace data. It also includes a generator for multivariate, 
# correlated count data. This generator will be used here to simulate data 
# based on existing data.
library("SpiecEasi")
set.seed(1234)
n <- 10 # number of samples

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

simulateNegBin <- function(phyloseqData, n, correlation=FALSE)
  
{ ##############################################################################
  # phyloseqData: 
  #   A phyloseq object containing the otu table which is to simulated
  #   using negative binomial distribution
  # n:
  #   The number of samples to be simulated
  # correlation:
  #   True if the count data needs to contain the correlation structure from the 
  #   input data
  #
  ##############################################################################
  # Output: count data drawn from a negative binomial distribution 
  #         n x p matrix with n the number of samples and p the number of otu's
  #
  ##############################################################################
  
  # Get otu table (DXN) and number of samples
  otuData <- otu_table(phyloseqData)
  
  # Transpose this to get rows to be individual samples
  X <- t(otuData)
  
  if (correlation)
  {
    # Simulate data using the otu means and otu variance covariance matrix from 
    # the input data, the shape parameter is estimated from the means and variance
    # matrix
    otuMeans <- colMeans(X)
    otuVarCovar <- var(X)
    simulatedData <- rmvnegbin(n = n, mu = otuMeans, Sigma = otuVarCovar)
  }
  else{
    # Size parameter per otu is estimated from the mean and variance
    # In case some are negative or zero, stop with an error message
    # Simulate data from single variate negative binomial distribution
    simulatedData <- apply(X, 2, rNegBin, n=n)
    row.names(simulatedData) <- 1:n
  }
  
  # Return the simulated data
  simulatedData
  
}

## THIS IS SOME CODE TO TEST THE FUNCTION, SHOULD BE REMOVED IN FINAL VERSION
testN <- 10
nrOtu <- 1000
nrSample <- 1000

otumat = matrix(sample(1:100, nrOtu*nrSample, replace = TRUE), nrow = nrOtu, ncol = nrSample)
rownames(otumat) <- paste0("OTU", 1:nrow(otumat))
colnames(otumat) <- paste0("Sample", 1:ncol(otumat))

taxmat = matrix(sample(letters, 70, replace = TRUE), nrow = nrow(otumat), ncol = 7)
rownames(taxmat) <- rownames(otumat)
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)

physeq = phyloseq(OTU, TAX)



testSim <- simulateNegBin(phyloseqData = physeq1, n=testN, correlation=TRUE)

###############################################################################


  








                      
                        