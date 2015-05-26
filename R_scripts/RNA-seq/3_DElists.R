# Assumes script 2_MergeDEtables.R was run previously


# General set up ----------------------------------------------------------


# Folder where the results of the analysis will be saved
rootDir = 'C:/Users/krue/Documents/Kevin-Logs/2013-08-27_AlvMac_RNAseq/04_complete_analysis'

setwd(rootDir)


# Create a list of DE features at each contrast ---------------------------


# Import the merged DE table produced previously
load('all.contrasts.ensembl71.RData')

# Initialise the super-list
DElists = list()

# Add the DE Ensembl gene identifiers
DElists$MB.CN.2 = all.contrasts.annotated$ensembl_gene_id[
  which(all.contrasts.annotated$FDR.MB.CN.2H < 0.05)]
DElists$MB.CN.6 = all.contrasts.annotated$ensembl_gene_id[
  which(all.contrasts.annotated$FDR.MB.CN.6H < 0.05)]
DElists$MB.CN.24 = all.contrasts.annotated$ensembl_gene_id[
  which(all.contrasts.annotated$FDR.MB.CN.24H < 0.05)]
DElists$MB.CN.48 = all.contrasts.annotated$ensembl_gene_id[
  which(all.contrasts.annotated$FDR.MB.CN.48H < 0.05)]

DElists$TB.CN.2 = all.contrasts.annotated$ensembl_gene_id[
  which(all.contrasts.annotated$FDR.TB.CN.2H < 0.05)]
DElists$TB.CN.6 = all.contrasts.annotated$ensembl_gene_id[
  which(all.contrasts.annotated$FDR.TB.CN.6H < 0.05)]
DElists$TB.CN.24 = all.contrasts.annotated$ensembl_gene_id[
  which(all.contrasts.annotated$FDR.TB.CN.24H < 0.05)]
DElists$TB.CN.48 = all.contrasts.annotated$ensembl_gene_id[
  which(all.contrasts.annotated$FDR.TB.CN.48H < 0.05)]

DElists$MB.TB.2 = all.contrasts.annotated$ensembl_gene_id[
  which(all.contrasts.annotated$FDR.MB.TB.2H < 0.05)]
DElists$MB.TB.6 = all.contrasts.annotated$ensembl_gene_id[
  which(all.contrasts.annotated$FDR.MB.TB.6H < 0.05)]
DElists$MB.TB.24 = all.contrasts.annotated$ensembl_gene_id[
  which(all.contrasts.annotated$FDR.MB.TB.24H < 0.05)]
DElists$MB.TB.48 = all.contrasts.annotated$ensembl_gene_id[
  which(all.contrasts.annotated$FDR.MB.TB.48H < 0.05)]

# Save the object to file
save(DElists, file = 'DElists.RData')
