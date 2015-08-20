
# Assumes script 9_GOexpress.R was run previously


# Libraries and dependencies ----------------------------------------------

library(GOexpress)


# General set up ----------------------------------------------------------


# Folder where the results of the analysis will be saved
rootDir = 'C:/Users/krue/Documents/Kevin-Logs/2013-08-27_AlvMac_RNAseq/04_complete_analysis'

setwd(rootDir)

# Create a folder to store raw counts
output.folder = 'ALOX-PTGS2'
dir.create(file.path(rootDir, output.folder))

# Load the ExpressionSet and GOexpress results
load('alvmac.eSet.RData')
load('GOexpress/gox.RData')

# Examine the phenotypic data
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


# Expression plot of the gene ALOX5 in various ways
# "smooth" plot restricted to 2-48 hpi (0 messes things, see below)
pdf(file = file.path(output.folder, 'ALOX5.plot_2-48.pdf'), width = 8, height = 6)
expression_plot_symbol(
  gene_symbol = 'ALOX5', result = gox, eSet = alvmac.eSet, x_var = 'Timepoint',
  subset=list(Time=c('2H', '6H', '24H', '48H')))
dev.off()
# "smooth" including 0 hpi which messes the smoothing
pdf(file = file.path(output.folder, 'ALOX5.plot_0-48.pdf'), width = 8, height = 6)
expression_plot_symbol(
  gene_symbol = 'ALOX5', result = gox, eSet = alvmac.eSet, x_var = 'Timepoint')
dev.off()
# Individual sample series
pdf(file = file.path(output.folder, 'ALOX5.plot_individual.pdf'), width = 8, height = 6)
expression_profiles_symbol(
  gene_symbol = 'ALOX5', result = gox, eSet = alvmac.eSet, x_var = 'Timepoint',
  seriesF = 'Series')
dev.off()


# Expression plot of the gene PTGS2 (BT.104982 in our data set)
# "smooth" plot restricted to 2-48 hpi (0 messes things, see below)
pdf(file = file.path(output.folder, 'PTGS2.plot_2-48.pdf'), width = 8, height = 6)
expression_plot_symbol(
  gene_symbol = 'BT.104982', result = gox, eSet = alvmac.eSet, x_var = 'Timepoint',
  subset=list(Time=c('2H', '6H', '24H', '48H')))
dev.off()
# "smooth" including 0 hpi which messes the smoothing
pdf(file = file.path(output.folder, 'PTGS2.plot_0-48.pdf'), width = 8, height = 6)
expression_plot_symbol(
  gene_symbol = 'BT.104982', result = gox, eSet = alvmac.eSet, x_var = 'Timepoint')
dev.off()
# Individual sample series
pdf(file = file.path(output.folder, 'PTGS2.plot_individual.pdf'), width = 8, height = 6)
expression_profiles_symbol(
  gene_symbol = 'BT.104982', result = gox, eSet = alvmac.eSet, x_var = 'Timepoint',
  seriesF = 'Series')
dev.off()
