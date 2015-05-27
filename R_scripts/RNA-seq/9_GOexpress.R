
# Assumes script 8_CreateExpressionSet.R was run previously


# Libraries and dependencies ----------------------------------------------


library(biomaRt)
library(ggplot2)
library(GOexpress)
library(tools)


# General set up ----------------------------------------------------------


# Folder where the results of the analysis will be saved
rootDir = 'C:/Users/krue/Documents/Kevin-Logs/2013-08-27_AlvMac_RNAseq/04_complete_analysis'

setwd(rootDir)

# Create a folder to store the volcano plots
goexpress.dir = 'GOexpress'
dir.create(file.path(rootDir, goexpress.dir))

# Load te 
load("alvmac.eSet.RData")


# Manually collect gene annotations from Ensembl release 71 ---------------


# Use ensembl mart with bovine dataset of Ensembl 71
ensembl71 <- useMart(
  host = 'apr2013.archive.ensembl.org',
  biomart='ENSEMBL_MART_ENSEMBL',
  dataset="btaurus_gene_ensembl"
  )

# Gene annotations
allgenes.ensembl71 <- getBM(
  attributes = c('ensembl_gene_id', 'external_gene_id', 'description'),
  mart = ensembl71)
colnames(allgenes.ensembl71)[1] <- 'gene_id'
head(allgenes.ensembl71)
# Save the object to file
save(
  allgenes.ensembl71,
  file=file.path(goexpress.dir, 'allgenes.ensembl71.RData')
  )

# Gene ontology annotations
allgo.ensembl71 <- getBM(
  attributes = c('go_id', 'name_1006', 'namespace_1003'),
  mart = ensembl71)
head(allgo.ensembl71)
# Save the object to file
save(
  allgo.ensembl71,
  file=file.path(goexpress.dir, 'allgo.ensembl71.rda')
  )

# Mapping of genes to gene ontologies
GOgenes.ensembl71 <- getBM(
  attributes = c('ensembl_gene_id', 'go_id'),
  mart = ensembl71)
colnames(GOgenes.ensembl71)[1] <- 'gene_id'
sum(GOgenes.ensembl71$gene_id == '')
sum(GOgenes.ensembl71$go_id == '')
GOgenes.ensembl71 <- GOgenes.ensembl71[GOgenes.ensembl71$go_id != '',]
head(GOgenes.ensembl71)
# Save the object to file
save(
  GOgenes.ensembl71,
  file=file.path(goexpress.dir, 'GOgenes.ensembl71.rda')
  )

# Remove the mart object that is no needed anymore
rm(ensembl71)


# GOexpress differentiating Infection -------------------------------------


# Let us look for genes that best classify the different infections groups
# across the time course of the experiment

# Run the GOexpress analysis
gox <- GO_analyse(
  eSet = alvmac.eSet, f = 'Infection',
  GO_genes = GOgenes.ensembl71,
  all_GO = allgo.ensembl71,
  all_genes = allgenes.ensembl71
  )

# Save the object to file
save(
  gox,
  file=file.path(goexpress.dir, 'gox.RData')
  )

# Further compress the file to optimise disk space usage
resaveRdaFiles(file.path(goexpress.dir, 'gox.RData'))

# Examine the gene and gene ontology ranking tables
head(subset_scores(gox, total=6)$GO, n = 20)
head(gox$genes, n=10)

# Remove the annotations that we don't need anymore
rm(GOgenes.ensembl71, allgenes.ensembl71, allgo.ensembl71)

# Export the gene ranking table
write.table(
  x = cbind(
    ensembl_gene_id=rownames(gox.pvalue$genes),
    gox.pvalue$genes),
  file = file.path(goexpress.dir, 'gox.pvalue.genes.txt'),
  quote = F, sep = '\t', row.names = F)


# GOexpress P-value -------------------------------------------------------

# Compute P-values for GO terms using random permutations of the genes
gox.pvalue <- pValue_GO(result = gox, N = 1000)

# Save the object to file
save(
  gox.pvalue,
  file=file.path(goexpress.dir, 'gox.pvalue.rda')
  )

# Further compress the file to optimise disk space usage
resaveRdaFiles(file.path(goexpress.dir, 'gox.pvalue.rda'))

# Export the GO term ranking table
write.table(
  x = gox.pvalue$GO,
  file = file.path(goexpress.dir, 'gox.pvalue.GO.txt'),
  sep = '\t', row.names = F, quote = F)


# Sample gene expression profiles -----------------------------------------


# Examine the top ranking gene list
head(gox$genes, n=10)

# Examine the phenotype information available to group samples
head(pData(alvmac.eSet))

# Create a numeric (i.e. not a R factor) representation of the time point
# information
# This will be used to plot the expression profile across the time course
# of the experiment, while respective the relative distance between the
# various time points.
alvmac.eSet$Timepoint = as.numeric(
  as.character(
    gsub('H', '', alvmac.eSet$Time)
    )
  )

# Expression plot of the top ranking gene SAA3 across all time points
# post-infection
pdf(file = file.path(goexpress.dir, 'SAA3.plot.pdf'), width = 8, height = 6)
expression_plot_symbol(
  gene_symbol = 'SAA3', result = gox, eSet = alvmac.eSet, x_var = 'Timepoint',
  subset=list(Time=c('2H', '6H', '24H', '48H')))
dev.off()

# Individual expression profiles of the top ranking gene SAA3 across all time points
# post-infection
pdf(file = file.path(goexpress.dir, 'SAA3.profiles.pdf'), width = 8, height = 6)
expression_profiles_symbol(
  gene_symbol = 'SAA3', result = gox, eSet = alvmac.eSet, x_var = 'Timepoint',
  seriesF = 'Series',
  subset=list(Time=c('2H', '6H', '24H', '48H')))
dev.off()


# Sample GO heatmaps ------------------------------------------------------


# Examine the top ranking list of biological processes
head(
  subset_scores(
    gox,
    total = 6,
    namespace = 'BP'
    )$GO,
  n = 20
  )

# Heatmap of the top GO term after filtering
# GO:0032727
# positive regulation of interferon-alpha production

# For clarity, plot one heatmap by time point
# (Increasing clustering of MB and CN, with TB as an intermediate
# phenotype)

pdf(
  file = file.path(goexpress.dir, 'IFNa_heatmap_2h.pdf'),
  width = 8, height = 6)
heatmap_GO(
  go_id = 'GO:0032727', result = gox, eSet = alvmac.eSet,
  subset = list(Time=c('2H')), cexCol = 1.5, cexRow = 0.9,
  main.Lsplit = 40
  )
dev.off()

pdf(
  file = file.path(goexpress.dir, 'IFNa_heatmap_6h.pdf'),
  width = 8, height = 6)
heatmap_GO(
  go_id = 'GO:0032727', result = gox, eSet = alvmac.eSet,
  subset = list(Time=c('6H')), cexCol = 1.5, cexRow = 0.9,
  main.Lsplit = 40
)
dev.off()

pdf(
  file = file.path(goexpress.dir, 'IFNa_heatmap_24h.pdf'),
  width = 8, height = 6)
heatmap_GO(
  go_id = 'GO:0032727', result = gox, eSet = alvmac.eSet,
  subset = list(Time=c('24H')), cexCol = 1.5, cexRow = 0.9,
  main.Lsplit = 40
)
dev.off()

pdf(
  file = file.path(goexpress.dir, 'IFNa_heatmap_48h.pdf'),
  width = 8, height = 6)
heatmap_GO(
  go_id = 'GO:0032727', result = gox, eSet = alvmac.eSet,
  subset = list(Time=c('48H')), cexCol = 1.5, cexRow = 0.9,
  main.Lsplit = 40
)
dev.off()


# Sample scoring table within GO term -------------------------------------

# Individual scores for the GO term
# GO:0032727
# positive regulation of interferon-alpha production
genes = table_genes(
  go_id = 'GO:0032727', result = gox.pvalue, data.only = F
  )

# Examine the individual gene scores
genes

# Export the table to a text file
write.table(
  x = cbind(
    ensembl_gene_id=rownames(genes),
    genes),
  file = file.path(goexpress.dir, 'IFNa_table.txt'),
  quote = F,
  sep = '\t', row.names = F)

# Remove the temporary object
rm(genes)


# Prepare data files for a Shiny application ------------------------------

# Create a sub-folder for the Shiny application
shiny.dir = 'shiny'
dir.create(file.path(goexpress.dir, shiny.dir))

# Create a sub-folder for the Shiny data files
data.dir = 'data'
dir.create(file.path(goexpress.dir, shiny.dir, data.dir))

# Rename the levels of the 'Infection' factor to human-friendly labels
alvmac.eSet$Infection = factor(
  alvmac.eSet$Infection, levels = c('CN', 'TB', 'MB'),
  labels = c('Control', 'M. tuberculosis', 'M. bovis'))

# Rename the levels of the 'Time' factor to human-friendly labels
alvmac.eSet$Time = factor(
  alvmac.eSet$Time, levels = c('0H', '2H', '6H', '24H', '48H'),
  labels = c('0 h', '2 hrs', '6 hrs', '24 hrs', '48 hrs'))

# Rename the 'Timepoint' column to 'Hours post-infection'
colnames(pData(alvmac.eSet))[colnames(pData(alvmac.eSet)) == 'Timepoint'] <- 'Hours post-infection'

# Save the edited expression set
saveRDS(
  alvmac.eSet,
  file = file.path(goexpress.dir, shiny.dir, data.dir, 'alvmac.eSet.rds')
  )

# List alphabetically non-empty gene names present in the expression set
# These will be used to populate the a dropdown menu in the application
genes_choices = sort(unique(gox$genes$external_gene_name))
genes_choices = genes_choices[genes_choices != '']
saveRDS(
  sort(unique(gox$genes$external_gene_name)),
  file = file.path(goexpress.dir, shiny.dir, data.dir, 'external_gene_names.rds')
  )

# Remove the temporary object
rm(genes_choices)

# List alphabetically non-empty GO terms present in the expression set
# These will be used to populate the a dropdown menu in the application
go_choices <- gox$GO$go_id
names(go_choices) <- gox$GO$name_1006
go_choices <- as.list(go_choices)
go_choices <- go_choices[sort(names(go_choices))]
saveRDS(
  go_choices,
  file = file.path(goexpress.dir, shiny.dir, data.dir, 'go_choices.rds')
  )

# Remove the temporary object
rm(go_choices)

# Save the GOexpress result object (including P-values)
saveRDS(
  gox.pvalue,
  file = file.path(goexpress.dir, shiny.dir, data.dir, 'gox.pvalue.rds')
  )

# Pre-extract the Ensembl gene identifiers present in the expression set
saveRDS(
  sort(rownames(alvmac.eSet)),
  file = file.path(goexpress.dir, shiny.dir, data.dir, 'ensemblIDs.rds')
)

# Pre-extract the time points present in the expression set
saveRDS(
  levels(alvmac.eSet$Time),
  file = file.path(goexpress.dir, shiny.dir, data.dir, 'timepoints.rds')
)

# Pre-extract the animal identifiers present in the expression set
saveRDS(
  sort(levels(alvmac.eSet$Animal)),
  file = file.path(goexpress.dir, shiny.dir, data.dir, 'animalIDs.rds')
)

# The above 'data/' folder and RDS files were subsequently copied
# into the 'R_scripts/shiny-alvmac' folder of this repository
# to complete the Shiny application


# Cleanup -----------------------------------------------------------------


# Remove all the objects saved to file or temporary in all the above sections
rm(alvmac.eSet, data.dir, goexpress.dir, gox, gox.pvalue, shiny.dir)


