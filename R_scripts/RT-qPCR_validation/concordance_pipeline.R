

# Libraries and dependencies ----------------------------------------------


# Load the required packages
library(ggplot2)
library(grid)
library(VennDiagram)
library(RColorBrewer)
library(magrittr)
library(gdata)
library(psych)
library(limma)


# General set up ----------------------------------------------------------

# Folder where the RT-qPCR data are stored
RTqPCR.folder = 'C:/Users/krue/Dropbox/Supervisors/Thesis/Chapter3_AlvMac/RT-qPCR_manual_Ct(20042015)'

# Folder where the results of the analysis will be saved
rootDir = 'C:/Users/krue/Documents/Kevin-Logs/2013-08-27_AlvMac_RNAseq/04_complete_analysis/RT-qPCR'
dir.create(rootDir)
setwd(rootDir)


# Import, filter, and annotate the RT-qPCR data ---------------------------


# Read in the excel RT-qPCR expression file
PCR <- read.xls(
  xls = file.path(RTqPCR.folder, 'CNRQ_tech_val.xlsx'),
  sheet = 1, row.names = 1, header = TRUE, na.strings = "NaN"
  )

# Clean up the column name in the qPCR dataset
colnames(PCR) %<>% gsub(
  pattern = "\\.SE.*", replacement = "_SE", x = .,
  perl = TRUE) %>% gsub(
    pattern = "\\.CNRQ.*", replacement = "_CNRQ", x = .,
    perl = TRUE) %>% gsub(
      pattern = "\\.", replacement = "_", x = .,
      perl = TRUE)

# Examine the edited table
PCR[1:5,1:5]

# Filter samples of interest
PCR.filtered = PCR[
  grep(pattern = 'W|NTC|STD|RT|\\.\\d$', x = rownames(PCR), invert = TRUE),]

# Remove the raw object
rm(PCR)

# Filter out the samples not done in technical validation
PCR.filtered = PCR.filtered[
  grep(
    pattern = 'N205|N206|N120|N178|N177', x = rownames(PCR.filtered),
    invert = TRUE),
  ]

# Extract phenotypic information
sample = data.frame(strsplit2(x = rownames(PCR.filtered), split = ' '))

# Format the phenotypic information
colnames(sample) <- c("animal", "infection", "timepoint")
rownames(sample) = paste(
  sample$animal, sample$infection, sample$timepoint, sep = ' ')

# Examine the edited table
head(sample)

# Merge the sample information with the PCR data
PCR.data <- merge(
  x = sample, y = PCR.filtered, by = "row.names", all.x = TRUE)

# Restore the sample name to the row names instead of a data column
rownames(PCR.data) <- PCR.data$Row.names
PCR.data$Row.names = NULL

# Remove the non-annotated PCR results
rm(PCR.filtered)


# Calculate the log fold-change in expression -----------------------------


# Log 2 transform the CNRQ value
PCR.data.log2 <- data.frame(
  PCR.data[, 1:(ncol(sample))],
  log2(x = PCR.data[, grep(pattern = "_CNRQ", x = colnames(PCR.data))]),
  row.names = row.names(PCR.data)
)

# Remove the raw PCR data and phenotypic annotations
rm(PCR.data, sample)

# Rename the column to indicate the data is log-transformed
colnames(PCR.data.log2) <- gsub(
  pattern = "_(CNRQ)", replacement = "_log\\1",x = colnames(PCR.data.log2))

# Define the variables required to compute fold-change in expression
targets.animal <- unique(PCR.data.log2$animal)
targets.time <- unique(PCR.data.log2$timepoint)
# Order the levels of factor 'time'
targets.time = factor(
  x = c('0H','2H','6H','24H','48H'), levels = c('0H','2H','6H','24H','48H'),
  ordered = TRUE)
# The 
targets.infection <- unique(PCR.data.log2$infection)

# Generate data for "animal" field in the fold change table
animal.column <- lapply(
  X = targets.animal,
  FUN = function(x) rep(
    x = x, length(targets.time)*length(targets.infection))) %>% unlist()

# Generate data for "timepoint" field in the fold change table
time.column = rep(x = lapply(
  X = targets.time,
  FUN = function(x) rep(
    x = x,
    length(targets.infection))) %>% unlist(), length(targets.animal)
)

# Generate data for "infection_contrast" column in the fold change table
contrast.column = rep(
  c('TB-CN', 'MB-CN', 'MB-TB'),
  length(targets.animal)*length(targets.time)
  )

full.rownames <- paste(animal.column, time.column, contrast.column, sep = "_")

# Initialise the data frame that will contains all the fold change values
logFC <- data.frame(row.names = full.rownames)
logFC <- cbind(
  logFC,
  strsplit2(x = rownames(logFC), split = '_')
)
colnames(logFC) <- c("animal", "timepoint", "contrast")
head(logFC)

# Remove the temporary objects
rm(full.rownames)

# Remove contrasts that cannot be computed
# (due to samples not processed/dropped)
logFC.available = logFC[grep(
  pattern = 'N1861_48H|N121_6H.*TB|N98_6H.*TB', x = rownames(logFC),
  perl = T, invert = TRUE
),]

# Remove the unfiltered logFC table
rm(logFC)

# for each gene
for (gene in colnames(PCR.data.log2)[4:ncol(PCR.data.log2)]){
  col.data <- c()
  # for each contrast defined in the logFC table
  for (contrast in rownames(logFC.available)){
    # get the different elements of the contrast
    contrast.data = strsplit2(contrast, split = '_')
    animal = contrast.data[1]
    time = contrast.data[2]
    # for time 0H insert artificial zeros (no data available for TB and MB)
    if(time == '0H'){
      col.data <- c(col.data, 0)
      next
    }
    infection.contrast = strsplit2(contrast.data[3], split = '-')
    infection.ref = infection.contrast[2]
    infection.target = infection.contrast[1]
    # get the data for the reference infection
    data.ref <- PCR.data.log2[
      PCR.data.log2$timepoint == time &
        PCR.data.log2$animal == animal &
        PCR.data.log2$infection == infection.ref
      , gene]
    # get the data for the target infection
    data.target <- PCR.data.log2[
      PCR.data.log2$timepoint == time &
        PCR.data.log2$animal == animal &
        PCR.data.log2$infection == infection.target
      , gene]
    # 
    col.data <- c(col.data, data.target-data.ref)
  }
  logFC.available[, gsub(pattern = "CNRQ", replacement = "FC", x = gene)] <- col.data
}
# Remove the temporary objects used in the loop
rm(
  gene, col.data, contrast, contrast.data, animal, time, infection.contrast,
  infection.ref, infection.target, data.ref, data.target
  )

# Examine the result table
head(logFC.available)


# Remove the log2-transforemed PCR data
rm(PCR.data.log2)


# Assess normal distribution of PCR data for each gene --------------------


# Create a copy of logFC without the 0H time point
logFC.infected = logFC.available[
  grep(pattern = '0H', x = rownames(logFC.available), invert = TRUE),]

# Plot the Q-Q plots for the overall data
qqplot <- apply(
  X = logFC.infected[, 4:ncol(logFC.infected)],
  MARGIN = 2,
  FUN = function(x) qqnorm(y = x, main = colnames(x)))

# Remove the temporary objects
rm(logFC.infected, qqplot)


# Compute significance values of fold-change in expression ----------------


# keep ignoring the time zero
# we will add them at the end for aesthetic

# Prepare data frame to include all fold-change and significance evaluation
gene.list <- colnames(logFC.available)[4:ncol(logFC.available)] %>% 
  lapply(
    X = .,
    FUN = function(x) rep(
      x = x, length(targets.time)*length(targets.infection))
  ) %>% unlist()

# get rid of the "logFC" text to define the row names
full.rownames <- paste(
  gene.list, time.column, contrast.column, sep = ".") %>% gsub(
    pattern = "_logFC.", replacement = ".", x = .)

# Define the name of the output table columns
full.colnames <- c(
  "gene", "time_point", "contrast", "shapiro", "n", "meanlogFC", "median",
  "sd", "se", "min", "max", "t.Pvalue", "w.Pvalue",
  "final.Pvalue")

# Initialise a matrix to contain the results
sig <- matrix(nrow = length(full.rownames), ncol = length(full.colnames))

# Convert the matrix to a data frame with named rows
sig <- data.frame(x = sig, row.names = full.rownames)

# Name the columns of the data frame
colnames(sig) <- full.colnames

# Insert the phenotypic information in the appropriate columns
sig[, c("gene", "time_point", "contrast")] <- strsplit2(
  x = row.names(sig), split =  ".", fixed = TRUE)


# Remove temporary objects
rm(
  animal.column, contrast.column, time.column, full.colnames, full.rownames,
  gene.list, targets.animal, targets.infection, targets.time
  )


# Compute the Pvalues and fold-change -------------------------------------


# for each contrast defined in the logFC table
for (contrast in rownames(sig)){
  shap.eval = list()
  t.eval = list()
  w.eval = list()
  # get the different elements of the contrast
  contrast.data = strsplit2(contrast, split = '.', fixed=TRUE)
  gene = contrast.data[1]
  time = contrast.data[2]
  infection.contrast = contrast.data[3]
  # get the data for the target group
  logFC.target <- logFC.available[
    logFC.available$timepoint == time  &
      logFC.available$contrast == infection.contrast,
    paste(gene, '_logFC', sep='')]
  # get general stats about target data
  stat.value <- describe(logFC.target)
  # test normality of target data
  if (time == "0H" | sum(!is.na(logFC.target)) < 3) {
    # Ignore time 0H, because no fold change to process
    shap.eval$p.value <- list("NaN")
    t.eval$p.value <- list("NaN")
    w.eval$p.value <- list("NaN")
  }
  else {
    shap.eval <- shapiro.test(x = logFC.target)
    # test differential expression (different from 0) using t-test
    # artificially create paired "0 fold change" samples for control logFC
    t.eval <- t.test(
      x = logFC.target, y = rep(0, length(logFC.target)),
      alternative = "two.sided", paired = TRUE,
      conf.level = 0.95, na.action = omit())
    # test differential expression (different from 0) using wilcoxon
    # artificially create paired "0 fold change" samples for control logFC
    w.eval <- wilcox.test(
      x = logFC.target, y = rep(0, length(logFC.target)),
      alternative = "two.sided", paired=TRUE, na.action = omit())
  }
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
  # Add result data in final table
  sig[contrast,
      c("shapiro", "n", "meanlogFC", "median", "sd", "se", "min",
        "max", "t.Pvalue", "w.Pvalue", "final.Pvalue")
      ] <- c(
        shap.eval$p.value, stat.value$n, stat.value$mean,
        stat.value$median, stat.value$sd, stat.value$se,
        stat.value$min, stat.value$max, t.eval$p.value,
        w.eval$p.value, final.pvalue
      )
}

# Remove the temporary objects used in the loop above
rm(
  shap.eval, t.eval, w.eval, contrast.data, gene, time, infection.contrast,
  logFC.target, stat.value, final.pvalue, contrast, logFC.available
  )

# Add a column with stars for significance
sig$sig.symbol = sapply(
  X = sig$final.Pvalue,
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

sig

save(sig, file='sig.rda')




# Import the corresponding RNAseq data ------------------------------------


# Import the RNA-seq DE table
RNAseq = read.xls(
  xls = file.path(RTqPCR.folder, 'RNAseq_10genes.xlsx'), sheet = 1,
  row.names = NULL, header = TRUE)

# Build a RNAseq table in the same format as the RT-qPCR
logFCs = c()
FDRs = c()
for (contrast in rownames(sig)){
  # get the contrast data from its name
  contrast.data = strsplit2(x = contrast, split = '.', fixed = TRUE)
  gene = contrast.data[1]
  time = contrast.data[2]
  infection.contrast = strsplit2(x = contrast.data[3], split = '-')
  infection.target = infection.contrast[1]
  infection.ref = infection.contrast[2]
  # fetch the RNAseq  data corresponding to the contrast
  if (time == '0H'){
    RNAseq.logFC = 0
    RNAseq.FDR = 'NaN'
  } else {
    RNAseq.logFC = RNAseq[
      RNAseq$external_gene_id == gene,
      grep(pattern = paste('logFC',infection.target, infection.ref,time, sep = '.*'), x = colnames(RNAseq), perl = TRUE)
      ]
    RNAseq.FDR = RNAseq[
      RNAseq$external_gene_id == gene,
      grep(pattern = paste('FDR',infection.target, infection.ref,time, sep = '.*'), x = colnames(RNAseq), perl = TRUE)
      ]
  }
  if (length(RNAseq.logFC) > 1){
    stop('RNAseq.logFC of length', length(RNAseq.FC))
  }
  if (length(RNAseq.FDR) > 1){
    stop('RNAseq.FDR of length', length(RNAseq.FDR))
  }
  logFCs = c(logFCs, RNAseq.logFC)
  FDRs = c(FDRs, RNAseq.FDR)
}
# Remove the temporary objects used in the loop
rm(
  contrast, contrast.data, gene, time, infection.contrast, infection.ref,
  infection.target, RNAseq.logFC, RNAseq.FDR
  )

# Assemble the results into a table 
sig.RNAseq = data.frame(
  RNAseq.logFC=logFCs, RNAseq.FDR=FDRs, row.names = rownames(sig)
  )

# Cleanup the last temporary objects
rm(logFCs, FDRs, RNAseq)

# Merge RT-qPCR and RNAseq DE tables
sig.merged = merge(x = sig, y = sig.RNAseq, by = 'row.names')

# Remove the two individual DE tables
rm(sig, sig.RNAseq)

# Add a column for RNAseq significance symbol
sig.merged$RNA.sig.symbol = sapply(
  X = sig.merged$RNAseq.FDR,
  FUN = function(x){
    x = as.numeric(as.character(x))
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
write.table(x = sig.merged, file = 'sig.merged.txt', sep = '\t', row.names = F)

# Remove the DE tests done using 7 animals or less
# as well as the time 0H (artificial samples meant for plotting)
# and also the reference genes (PPIA and H3F3A)
sig.merged.8more = sig.merged[
  which(
    sig.merged$n >= 8 &
      sig.merged$time_point != '0H' &
      !sig.merged$gene %in% c('PPIA', 'H3F3A')),]

# Save the object to file
save(sig.merged.8more, file='sig.merged.8more.rda')

# Count the number of contrasts left for each gene
table(sig.merged.8more$gene)

# Note that we lost PIK3IP1 because all non-zero contrasts have 6 animals or less
# TLR2 and IL6 have data for less than the expected 12 contrasts per gene

# Concordance 
# Both tests are not significant
# Both test are significant with the same direction of fold-change
# 0.6621622 (66.2%)
sum(
  apply(
    X = sig.merged.8more[,c('meanlogFC','sig.symbol','RNAseq.logFC','RNA.sig.symbol')],
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
  nrow(sig.merged.8more)


# Correlation (all) Pearson 0.98
plot(x = sig.merged.8more$meanlogFC, y = sig.merged.8more$RNAseq.logFC)
cor.test(
  x = sig.merged.8more$meanlogFC,
  y = sig.merged.8more$RNAseq.logFC, method = 'pearson')
# cor = 0.9834018 - p < 2.2e-16
cor.test(
  x = sig.merged.8more$meanlogFC,
  y = sig.merged.8more$RNAseq.logFC, method = 'spearman')
# r = 0.9717734 - p < 2.2e-16

# Correlation (per gene)
cors = c()
Ps = c()
for (gene in unique(sig.merged.8more$gene)){
  cor_result = cor.test(
    x = sig.merged.8more[sig.merged.8more$gene == gene, 'meanlogFC'], # PCR
    y = sig.merged.8more[sig.merged.8more$gene == gene, 'RNAseq.logFC'], # RNAseq
    method = 'pearson')
  cors = c(cors, cor_result$estimate)
  Ps = c(Ps, cor_result$p.value)
}
# Remove the temporary objects used in the loop
rm(gene, cor_result)

genes.cor = data.frame(
  cor = cors,
  p.value = Ps,
  row.names = unique(sig.merged.8more$gene)
)
# Remove the objects now saved into genes.cor
rm(cors, Ps)
genes.cor

#            cor      p.value
# CCL4 0.9885694 1.507646e-09
# FOS  0.9707618 1.602211e-07
# IL10 0.9510528 2.037636e-06
# IL1B 0.9881214 1.825799e-09
# IL6  0.9737446 2.014189e-06
# TLR2 0.9761297 2.387031e-02
# TNF  0.9913064 3.854304e-10

# Write the individual correlation coefficients to a text file
write.table(
  x = genes.cor, file = 'genes.correlations.txt', sep = '\t', row.names = F,
  quote = FALSE)

# Remove the temporary and saved objects
rm(genes.cor, sig.merged)


# Concordance considering only MB_CN contrasts ----------------------------

# Checking MB-CN contrasts only
# (comparison with previous results obtained by Nicolas)

sig.MB_CN = sig.merged.8more[sig.merged.8more$contrast == 'MB-CN',]

sum(
  apply(
    X = sig.MB_CN[,c('meanlogFC','sig.symbol','RNAseq.logFC','RNA.sig.symbol')],
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
  nrow(sig.MB_CN)

# Remove the temporary object
rm(sig.MB_CN)



# Cleanup -----------------------------------------------------------------

# Remove the remaining temporary objects
rm(RTqPCR.folder, sig.merged.8more)
