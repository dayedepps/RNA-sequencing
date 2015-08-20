
# Libraries and dependencies ----------------------------------------------


library(Biobase)
library(gdata)
library(limma)
library(psych)


# General set up ----------------------------------------------------------


# Folder where the RT-qPCR data are stored
MSD.folder = 'C:/Users/krue/Dropbox/Supervisors/Thesis/Chapter3_AlvMac/MSD_protein(22042015)/'

# Folder where the results of the analysis will be saved
rootDir = 'C:/Users/krue/Documents/Kevin-Logs/2013-08-27_AlvMac_RNAseq/04_complete_analysis/MSD'
dir.create(rootDir)
setwd(rootDir)


# Load in the data --------------------------------------------------------


# Read in the excel RT-qPCR expression file
MSD <- read.xls(
  xls = file.path(MSD.folder, "MSD_raw_data(20042015)-KR.xlsx"),
  sheet = 2, row.names = 1, header = TRUE
)

# Examine the data
t(MSD)
dimnames(MSD)
dimnames(MSD)[c(2,1)]
matrix(t(MSD), nrow = ncol(MSD), ncol = nrow(MSD), dimnames = dimnames(MSD)[c(2,1)])
t(matrix(MSD, ))


# Create an ExpressionSet of the MSD data ---------------------------------


MSD.eset = ExpressionSet(
  assayData = matrix(
    t(MSD), nrow = ncol(MSD), ncol = nrow(MSD),
    dimnames = dimnames(MSD)[c(2,1)])
)
MSD.eset

# Annotate the ExpressionSet with sample phenodata
phenodata = data.frame(
  row.names = colnames(MSD.eset),
  strsplit2(x = colnames(MSD.eset), split = '_')
)
colnames(phenodata) = c('Animal', 'Infection', 'Time')
head(phenodata)

phenoData(MSD.eset) = AnnotatedDataFrame(phenodata)
head(pData(MSD.eset))

# Cleanup
rm(MSD, phenodata)


# Retain only the 10 animals processed in RNAseq ------------------------

# Restrict the data to the 10 animals processed in RNAseq
MSD.eset.10 = MSD.eset[
  ,grep(
    'N1178|N121|N138|N158R|N1855|N1859|N1861|N1864|N1870|N98',
    MSD.eset$Animal)
  ]

# Check the number of samples per animal
table(MSD.eset.10$Animal)
# For the animal that has no data at 48H, only 9 contrasts instead of 12
# We can also see the 5 extra animals not processed in RNAseq (0 samples)

MSD.eset.10

# Cleanup
rm(MSD.eset)


# List all contrasts to perform -------------------------------------------


# Re-order for convenience
MSD.eset.10$Time = factor(
  MSD.eset.10$Time,
  levels = c('2H', '6H', '24H', '48H'),
  ordered = TRUE
  )
levels(MSD.eset.10$Time)

# For each protein, 4 timepoints, and
# For each time-point, three contrasts
protein.column = unlist(
  sapply(
    X = rownames(MSD.eset.10),
    FUN = function(x) rep(x, length(levels(MSD.eset.10$Time))*3),
    simplify = F, USE.NAMES = FALSE)
)

time.column = rep(
  unlist(
    sapply(
      X = levels(MSD.eset.10$Time),
      FUN = function(x) rep(x, 3),
      simplify = F, USE.NAMES = FALSE)
  ),
  length(rownames(MSD.eset.10)))

contrast.column = rep(
  c('TB-CN', 'MB-CN', 'MB-TB'),
  length(levels(MSD.eset.10$Time))*length(rownames(MSD.eset.10))
)

full.rownames = paste(protein.column, time.column, contrast.column, sep='_')

# Assemble the result table
sig.table = data.frame(
  Protein=protein.column, Time=time.column, Contrast=contrast.column,
  row.names = full.rownames)

# Cleanup
rm(protein.column, time.column, contrast.column, full.rownames)


# Collect the statistics for all contrasts --------------------------------

# Initialise a temporary container for the results
stats = list(
  shap.eval.p.value=c(), stat.value.n=c(), stat.value.mean=c(),
  stat.value.median=c(), stat.value.sd=c(), stat.value.se=c(),
  stat.value.min=c(), stat.value.max=c(), t.eval.p.value=c(),
  w.eval.p.value=c(), final.pvalue=c())

# for each contrast defined in the logFC table
for (contrast in rownames(sig.table)){
  # get the different elements of the contrast
  contrast.data = strsplit2(contrast, split = '_', fixed=TRUE)
  protein = contrast.data[1]
  time = contrast.data[2]
  infection.contrast = strsplit2(contrast.data[3], split = '-', fixed=TRUE)
  infecton.target = infection.contrast[1]
  infecton.ref = infection.contrast[2]
  # Get the reference data
  ref.data = exprs(
    MSD.eset.10[
      protein,
      which(
        MSD.eset.10$Time == time &
          MSD.eset.10$Infection == infecton.ref
      )
      ]
  )
  # Get the target data
  target.data = exprs(
    MSD.eset.10[
      protein,
      which(
        MSD.eset.10$Time == time &
          MSD.eset.10$Infection == infecton.target
      )
      ]
  )
  # Residuals (paired)
  residuals = c(target.data - ref.data)
  # get general stats about residuals
  stat.value <- describe(residuals)
  # Test normality of the residuals
  shap.eval <- shapiro.test(residuals)
  #
  t.eval <- t.test(
    x = target.data, y = ref.data, 
    alternative = "two.sided", paired = TRUE,
    conf.level = 0.95)
  #
  w.eval <- wilcox.test(
    x = target.data, y = ref.data, 
    alternative = "two.sided", paired=TRUE)
  # Define p-value to use depending on the shapiro test
  if (shap.eval$p.value == "NaN") {
    final.pvalue <- "NaN"      
  }
  else if (shap.eval$p.value < 0.1) {
    final.pvalue <- w.eval$p.value
  }
  else {
    final.pvalue <- t.eval$p.value
  }
  stats[['shap.eval.p.value']] = c(stats[['shap.eval.p.value']], shap.eval$p.value)
  stats[['stat.value.n']] = c(stats[['stat.value.n']], stat.value$n)
  stats[['stat.value.mean']] = c(stats[['stat.value.mean']], stat.value$mean)
  stats[['stat.value.median']] = c(stats[['stat.value.median']], stat.value$median)
  stats[['stat.value.sd']] = c(stats[['stat.value.sd']], stat.value$sd)
  stats[['stat.value.se']] = c(stats[['stat.value.se']], stat.value$se)
  stats[['stat.value.min']] = c(stats[['stat.value.min']], stat.value$min)
  stats[['stat.value.max']] = c(stats[['stat.value.max']], stat.value$max)
  stats[['t.eval.p.value']] = c(stats[['t.eval.p.value']], t.eval$p.value)
  stats[['w.eval.p.value']] = c(stats[['w.eval.p.value']], w.eval$p.value)
  stats[['final.pvalue']] = c(stats[['final.pvalue']], final.pvalue)
}
# Remove the temporary variables used in the loop
rm(
  contrast, contrast.data, protein, time, infection.contrast,
  infecton.ref, infecton.target, ref.data, target.data, residuals,
  stat.value, shap.eval, t.eval, w.eval, final.pvalue)

# Examine the results
str(stats)
# 84 values in each vector

# Fill the corresponding columns of the result table
sig.table$shapiro = stats[['shap.eval.p.value']]
sig.table$n = stats[['stat.value.n']]
sig.table$meanlogFC = stats[['stat.value.mean']]
sig.table$median = stats[['stat.value.median']]
sig.table$sd = stats[['stat.value.sd']]
sig.table$se = stats[['stat.value.se']]
sig.table$min = stats[['stat.value.min']]
sig.table$max = stats[['stat.value.max']]
sig.table$t.Pvalue = stats[['t.eval.p.value']]
sig.table$w.Pvalue = stats[['w.eval.p.value']]
sig.table$final.Pvalue = stats[['final.pvalue']]

# Add a column with stars for significance
sig.table$sig.symbol = sapply(
  X = sig.table$final.Pvalue,
  FUN = function(x){
    if (x < 0.001){
      return('***')
    } else if (x < 0.01){
      return('**')
    } else if (x < 0.05){
      return('*')
    } else{
      return('')
    }
  })

save(sig.table, file='sig.table.rda')

# Cleanup
rm(stats, MSD.eset.10)


# Import the RNAseq data for the associated genes -------------------------


# Import the pre-filtered RNAseq data
RNAseq = read.xls(
  xls = file.path(MSD.folder, "RNAseq_6genes(no-BSA).xlsx"),
  sheet = 1, row.names = 1, header = TRUE
)

# # Check why ALB (gene for BSA) is absent from RNAseq
# load('C:/Users/krue/Documents/Kevin-Logs/2013-08-27_AlvMac_RNAseq/01_SenseCountsAnalysis/raw.RData')
# range(raw['ENSBTAG00000017121',]$counts)
# # only 0-12 raw reads per library

# Build a RNAseq table in the same format as the RT-qPCR
logFCs = c()
FDRs = c()
for (contrast in rownames(sig.table)){
  # get the contrast data from its name
  contrast.data = strsplit2(x = contrast, split = '_')
  gene = toupper(gsub('\\.', '', contrast.data[1])) # protein to gene conversion
  # Skip protein 'BSA'
  if (gene == 'BSA'){ # ALB gene was filtered out for low expression (0-12 reads)
    logFCs = c(logFCs, 'NaN')
    FDRs = c(FDRs, 'NaN')
    next
  } else if (gene == 'IL12'){ # protein to gene conversion
    gene = 'IL12B'
  } else if (gene == 'TNFA'){ # protein to gene conversion
    gene = 'TNF'
  }
  time = contrast.data[2]
  infection.contrast = strsplit2(x = contrast.data[3], split = '-')
  infection.target = infection.contrast[1]
  infection.ref = infection.contrast[2]
  # fetch the RNAseq  data corresponding to the contrast
  RNAseq.logFC = RNAseq[
    RNAseq$external_gene_id == gene,
    grep(
      pattern = paste('logFC',infection.target, infection.ref,time, sep = '.*'),
      x = colnames(RNAseq), perl = TRUE)
    ]
  RNAseq.FDR = RNAseq[
    RNAseq$external_gene_id == gene,
    grep(
      pattern = paste('FDR',infection.target, infection.ref,time, sep = '.*'),
      x = colnames(RNAseq), perl = TRUE)
    ]
  if (length(RNAseq.logFC) != 1){
    stop('RNAseq.logFC of length: ', length(RNAseq.logFC))
  }
  if (length(RNAseq.FDR) != 1){
    stop('RNAseq.FDR of length: ', length(RNAseq.FDR))
  }
  logFCs = c(logFCs, RNAseq.logFC)
  FDRs = c(FDRs, RNAseq.FDR)
}
# Remove teh temporary objects used in the loop above
rm(
  RNAseq.FDR, RNAseq.logFC, infection.contrast, infection.ref,
  infection.target, time, gene, contrast.data, contrast)

# Assemble the results in a data frame
sig.RNAseq = data.frame(
  RNAseq.logFC=logFCs, RNAseq.FDR=FDRs, row.names = rownames(sig.table))

# Cleanup
rm(logFCs, FDRs)

# Merge RT-qPCR and RNAseq DE tables
sig.merged = merge(x = sig.table, y = sig.RNAseq, by = 'row.names')

# Add a column with stars for significance (RNAseq)
sig.merged$RNA.sig.symbol = sapply(
  X = as.numeric(as.character(sig.merged$RNAseq.FDR)),
  FUN = function(x){
    if (x == 'NaN'){
      return('')
    }
    else if (x < 0.001){
      return('***')
    } else if (x < 0.01){
      return('**')
    } else if (x < 0.05){
      return('*')
    } else{
      return('')
    }
  })

# Save the object to file
save(sig.merged, file='sig.merged.rda')

# Write the table to a text file
write.table(
  x = sig.merged, file = 'sig.merged.txt', sep = '\t', row.names = F,
  quote = FALSE)

# Remove the individual DE tables
rm(sig.table, sig.RNAseq)


# Concordance ignoring BSA (ALB) ------------------------------------------

sig.merged.noBSA = sig.merged[sig.merged$Protein != 'BSA',]
sig.merged.noBSA$RNAseq.logFC = as.numeric(as.character(sig.merged.noBSA$RNAseq.logFC))

# Concordance:
# Both tests are not significant
# Both test are significant with the same direction of fold-change

# 0.8055556 (80.6 %)
sum(
  apply(
    X = sig.merged.noBSA[,c('meanlogFC','sig.symbol','RNAseq.logFC','RNA.sig.symbol')],
    MARGIN = 1,
    FUN = function(x){
      # PCR is not significant
      if (x[[2]] == ''){
        # Both technologies do not show significant DE
        if (x[[4]] == ''){
          return(TRUE)
        }
        # PCR is not significant, while RNA-seq is significant
        else{
          return(FALSE)
        }
      }
      # PCR is significant
      else{
        # PCR is significant, while RNA-seq is not significant
        if (x[[4]] == ''){
          return(FALSE)
        }
        # Both are significant
        else{
          PCR = as.numeric(as.character(x[[1]]))
          RNAseq = as.numeric(as.character(x[[3]]))
          # Both technologies show the same direction of fold change
          # then the product of fold-change values is positive
          if ((PCR * RNAseq) > 0){
            return(TRUE)
          }
          else {
            return(FALSE)
          }
        }
      }
    }
    )
  ) / 
  nrow(sig.merged.noBSA)

# Saving tables to temporary objects allows rapid visualsisation in RStudio
tmp = sig.merged.noBSA[sig.merged.noBSA$final.Pvalue < 0.05,]
