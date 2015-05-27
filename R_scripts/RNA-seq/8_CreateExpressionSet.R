
# Assumes script 1_alv-mac-edgeR-pipeline.R was run previously


# Libraries and dependencies ----------------------------------------------


library(tools)


# General set up ----------------------------------------------------------


# Folder where the results of the analysis will be saved
rootDir = 'C:/Users/krue/Documents/Kevin-Logs/2013-08-27_AlvMac_RNAseq/04_complete_analysis'

setwd(rootDir)

# Load the expression data and phenotypic information
load('log2_cpm_norm.RData')
load('targets.RData')



# Create an ExpressionSet of the Biobase package --------------------------

# Store the two objects in a dedicated Bioconductor container
# This step performs a few controls and ensures that the experimental
# data and phenotypic information remain consistent.
# It also represents a well-defined data container used by several
# Bioconductor packages to ease the handling of data.
alvmac.eSet = ExpressionSet(assayData = log2_cpm_norm, phenoData = targets)

# Save the object to file
save(alvmac.eSet, file='alvmac.eSet.RData')

# Further compress the file to optimise disk space usage
resaveRdaFiles('alvmac.eSet.RData')

# Remove the objects now saved in the ExpressionSet
rm(log2_cpm_norm, targets)
