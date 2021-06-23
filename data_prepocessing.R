# This script handles all the preprocessing of the phyloseq data object

# Look only at OTU's with a library size larger than minLibSize
phyloseqData <- prune_taxa(taxa_sums(phyloseqData) >= minLibSize, phyloseqData)

preprocData <- function(phyloData, minLibSize)
{
  ##############################################################################
  # phyloData:
  #   A phyloseq object containing the real data.
  # minLibSize:
  #   The minimum library size for an otu to be kept.
  ##############################################################################
  # Output: 
  #   A preprocessed version of the input phyloseq object. 
  ##############################################################################
  
  # Look only at OTU's with a library size larger than minLibSize
  phyloData <- prune_taxa(taxa_sums(phyloData) >= minLibSize, phyloData)
  
  phyloData
}