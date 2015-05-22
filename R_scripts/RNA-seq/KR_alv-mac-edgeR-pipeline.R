
# Libraries ---------------------------------------------------------------

library(biomaRt)
library(Biobase)
library(edgeR)
library(ggplot2)
library(PerformanceAnalytics)
library(tools)


# General set up ----------------------------------------------------------

# Folder where the count files are stored
countsFolder = 'C:/Users/krue/Documents/Kevin-Logs/2013-08-27_AlvMac_RNAseq/Counts'

# Folder where the results of the analysis will be saved
rootDir = 'C:/Users/krue/Documents/Kevin-Logs/2013-08-27_AlvMac_RNAseq/04_complete_analysis'

setwd(rootDir)


# Reading raw data files into edgeR ---------------------------------------


# List the counts files identified through a filename pattern
# Note: All RNAseq files in the folder end in 'H' for 'hours'
list.files(
  path = countsFolder, pattern='N*H$', all.files = FALSE, full.names = FALSE,
  recursive = FALSE, ignore.case = FALSE)

# Store the list of counts files
files = dir(path=countsFolder, pattern='N*H$')
files

# Create a list defining which treatment:time each sample belongs to
groups = factor(
  sapply(
    X = files,
    FUN = function(x){
      paste(unlist(strsplit(x[1], '_'))[2:3], collapse = '_')
      }
    ),
  levels = c(
    'CN_0H','CN_2H','CN_6H','CN_24H','CN_48H','MB_2H','MB_6H','MB_24H',
    'MB_48H','TB_2H','TB_6H','TB_24H','TB_48H')
  )
groups

# Read and merge the set of count files
# Note: need to specify which columns in the featureCounts files that contain
# the tag names (i.e. Ensembl gene ID) and the counts/library size (library
# size contains info on all gener counts).
# These are in columns 1 and 3, respectively.
raw = readDGE(files=files, path=countsFolder, columns=c(1,3), group=groups)

# Save the R object to file
save(raw, file='raw.RData')

# Further compress the file to optimise disk space usage
resaveRdaFiles('raw.RData')

# Remove the 'files' and 'groups' object which are now saved in the 'raw' data
rm(files, groups)


# Create an annotated data frame to describe the data ---------------------


# Isolates the different levels of the different factors to set up an
# annotation data frame variable
animal = factor(sapply(strsplit(raw$samples$files, '_'), '[[', 1))
infection = factor(
  sapply(strsplit(raw$samples$files, '_'), '[[', 2),
  levels = c('CN', 'TB', 'MB'))
time = factor(
  sapply(strsplit(raw$samples$files, '_'), '[[', 3),
  levels = c('0H','2H','6H','24H','48H'),
  ordered=TRUE)
sample = gsub(pattern='_S', replacement='', x=raw$samples$files)

# Build a variable containing all information for the samples in the study
targets = AnnotatedDataFrame(
  data=data.frame(
    File=raw$samples$files,
    Sample=sample,
    Animal=animal,
    Infection=infection,
    Time=time,
    Group=raw$samples$group,
    Series=paste(animal, infection, sep = '_'),
    stringsAsFactors=F)
  )

# Examine the object
head(pData(targets))

# Save the R object to file
save(targets, file='targets.RData')

# Removes all objects that are now saved in targets
rm(sample, animal, infection, time)

# Examine the 'targets' variable
unique(targets$Animal)
length(unique(targets$Animal))
unique(targets$Group)
length(unique(targets$Group))
unique(targets$Infection)
unique(targets$Time)

# Examine the number of animals (i.e. biologicla replicates) in each group
# We know that animal ID N1861 was not sequenced at 48 hpi
aggregate(Animal~Group, data=pData(targets), FUN=unique)
# The three treatments at 48H have one animal less than the other groups.
# It is always the 7th animal
unique(targets$Animal)[[7]]
# which is indeed animal ID N1861


# examine the library sizes for all files in 'raw' DGElist ----------------

# Examine the DGEList object
head(raw$samples)
raw$counts[1:5,1:5]

pdf(file = 'raw.library.size.pdf', width = 10, height = 7)
# Plot the raw library size
plot(
  raw$samples$lib.size,
  main='Raw library size',
  xlab='Sample index', ylab=' Library size',
  ylim = c(0, max(raw$samples$lib.size)))
# Save the figure to file
dev.off()

# One of the libraries has unusually large library size

# identify the index of the outlier having a larger library size
which(raw$samples$lib.size == max(raw$samples$lib.size))
# Indicates sample 117
raw$samples[117,]
# Indicates sample N98_CN_2H, with a library size of 27,930,001 read pairs


# QC analysis of each individual raw library ------------------------------


# Obtain log2-transformed counts per million (CPM) values for each sample
# The log2 transformation allows to scale the data and subsequent plots
log2_cpm_raw = cpm(x=raw, log=TRUE, normalized.lib.sizes=TRUE)
log2_cpm_raw[1:5,1:5]

# Save the R object to file
save(log2_cpm_raw, file='log2_cpm_raw.RData')

# Further compress the file to optimise disk space usage
resaveRdaFiles('log2_cpm_raw.RData')

# Plot the distribution of raw log2 CPM per library
# Note: highlight the outlier sample 117 (identified above) in red
pdf(file = 'log2.cpm.raw.pdf')
plot(
  density(log2_cpm_raw[,1]),
  main='Density plot of log2(cpm) raw count per gene',
  lty=1, ylim=c(0,0.31),
  xlab='Log2(cpm) of count per gene', ylab='Density',
  col='black')
for (i in 2:127) {
  if(i != 117){lines(density(log2_cpm_raw[,i]), lty=1, col='black')}
}
lines(density(log2_cpm_raw[,117]), lty=1, lwd=1.5, col='red')
rm(i)
dev.off()

# Remove the temporary object containing the raw log2(cpm)
rm(log2_cpm_raw)

# It seems the larger library size increases the detection of lowly expressed
# genes, but does not visibly affect highly expressed genes


# Compare the correlation of the outlier raw reads with other 2H samples -----

load('log2_cpm_raw.RData')

# Subset the raw data to the CN_2H samples
log2.cpm.raw.2H = log2_cpm_raw[,grep(pattern='CN_2H', x=pData(targets)$Group)]

# Rename the columns with only the animal name (time and treatment are
# the same for all)
colnames(log2.cpm.raw.2H)
colnames(log2.cpm.raw.2H) = sapply(
  strsplit(colnames(log2.cpm.raw.2H), split="_"),
  "[[",
  1)

# Pairwise correlations of samples in group CN_2H
pdf(file = 'chart_correlation.pdf', width = 10, height = 10)
chart.Correlation(R=log2.cpm.raw.2H, method="pearson")
dev.off()

# The outlier library shows a barely reduced r-squared with the other
# libraries (0.96-0.98) relative to r-squared values between all other
# libraries (0.97-0.98).

# Remove the temporary objects
rm(log2.cpm.raw.2H, log2_cpm_raw)


# Filter low expressed and rRNA genes ------------------------------------


# The cpm() function in edgeR Returns counts per million from a DGEList
# or matrix object by dividing raw counts by library size (which can be
# normalized) and multiplying by one million.
cpm(raw)[1:5,1:5]


# Filter low expression tags.
# use tags that achieve one count per million for at least 10 libraries.
# 10 libraries is based on the number of biological replicates
filtered <- raw[rowSums(cpm(raw) > 1) >= 10,]

# Save the filtered counts object to file
save(filtered, file='filtered.RData')

# Further compress the file to optimise disk space usage
resaveRdaFiles('filtered.RData')

# Examine the filtered counts object
dim(filtered)
# 12.121 genes remain after filtering in the 127 samples

# Read in list of rRNA genes into edgeR
# This list was obtained by extracting Ensembl gene identifiers of the
# Bos_taurus.UMD3.1.73.gtf file
# (ftp://ftp.ensembl.org/pub/release-73/gtf/bos_taurus/)
# annotated to 'rRNA' in the feature column
rRNA = scan(
  'C:/Users/krue/Documents/Kevin-Logs/2013-08-27_AlvMac_RNAseq/rRNA/rRNA_ensembl.txt',
  what = character())

length(rRNA)
## There are 405 nuclear-encoded rRNA genes in the list
# Note that the two mtDNA rRNA genes have been retained in the data set

rRNA.filter = rownames(filtered$counts) %in% rRNA ## check how many rRNA genes are in the 'filtered' DGElist
summary(rRNA.filter)
# 404 rRNA genes have been filtered out using the low level expression
# filtering criteria above while 1 rRNA gene is left

# Filter out the last rRNA gene from the 'filtered' DGElist
filt_rm_rRNA <- filtered[!rRNA.filter,] 

# Examine the resulting DGEList
filt_rm_rRNA[1:5,1:5]
dim(filt_rm_rRNA)

# There are now 12,121 genes (rows) measured in 127 samples (columns)

# Remove the temporary objects
rm(filtered, rRNA.filter, rRNA)

## Re-compute the library sizes after filtering (as per edgeR users manual)
filt_rm_rRNA$samples$lib.size <- colSums(filt_rm_rRNA$counts)

# Save the updated object
save(filt_rm_rRNA, file="filt_rm_rRNA.RData")

# Further compress the file to optimise disk space usage
resaveRdaFiles('filt_rm_rRNA.RData')


# Compare the library size before and after filtering ---------------------

# library sizes before filtering
raw$samples$lib.size
# library sizes after filtering
filt_rm_rRNA$samples$lib.size
# Note: the numbers are different (but not massively different)

# Number of reads in each sample before (black) and after (red) filtering
pdf(file = 'raw.and.filtered.library.size.pdf', width = 10, height = 7)
plot(
  raw$samples$lib.size,
  pch = 16, main = "Raw and filtered libraries size",
  xlab = "Sample index", ylab="Library size",
  ylim = c(0, max(raw$samples$lib.size)))
points(filt_rm_rRNA$samples$lib.size, col="red")
dev.off()

# Number of reads removed by sample
pdf(file = 'proportion.filtered.library.size.pdf', width = 10, height = 7)
plot(
  (raw$samples$lib.size - filt_rm_rRNA$samples$lib.size) / raw$samples$lib.size*100,
  xlab="Sample index", ylab="Proportion of reads filtered out (%)",
  main="Propotion of libraries filtered",
  ylim = c(0, 1))
dev.off()

# Remove the raw counts object not used anymore below
rm(raw)


# Normalisation of data ---------------------------------------------------


# Normalization of data using trimmed mean of M-values (normalization based
# on RNA composition between libraries)
# Calculate normalisation factors for the filtered DGElist (filt_rm_rRNA). 
# With edgeR, counts are not transformed in any way after normalization,
# instead normalization will modify library size)
filt_rm_rRNA_norm <- calcNormFactors(filt_rm_rRNA)

# Save the object containing the normalisation factors
save(filt_rm_rRNA_norm, file="filt_rm_rRNA_norm.RData")

# Further compress the file to optimise disk space usage
resaveRdaFiles('filt_rm_rRNA_norm.RData')

# Remove the pre-normalisation DGEList
rm(filt_rm_rRNA)


# QC post normalisation -------------------------------------------


# CPM values considering the normalisation factor for each library
log2_cpm_norm = cpm(x=filt_rm_rRNA_norm, log=TRUE, normalized.lib.sizes=TRUE)

# Save the R object to file
save(log2_cpm_norm, file="log2_cpm_norm.RData")

# Further compress the file to optimise disk space usage
resaveRdaFiles('log2_cpm_norm.RData')

# Plot the distribution of normalised log2 CPM per library
# Note: highlight the outlier sample 117 (identified above) in red
pdf(file = 'log2.cpm.normalised.pdf')
plot(
  density(log2_cpm_norm[,1]),
  main="Density plot of log2(cpm) normalised count per gene",
  ylim=c(0,0.18), lty=1,
  xlab="Log2(cpm) of count per gene", ylab="Density", col="black")
for (i in 2:127) {
  if(i != 117){lines(density(log2_cpm_norm[,i]), lty=1, col="black")}
}
lines(density(log2_cpm_norm[,117]), lty=1, lwd=1.5, col="red")
dev.off()
rm(i)

# Again, it seems the larger library size increases the detection of lowly
# expressed genes, but does not visibly affect highly expressed genes

# Remove the temporary object containing the normalised log2(cpm)
rm(log2_cpm_norm)


# MDS plot (general) ------------------------------------------------------

# Save the MDS object obtained using the limma package 
MDS = plotMDS(filt_rm_rRNA_norm, top=nrow(filt_rm_rRNA_norm))

# Save the object to file
save(MDS, file="MDS.RData")

# Extract the coordinates for each sample
MDS.data = data.frame(MDS$cmdscale.out)

# Check that samples are in the same order than the 'targets' 
# annotated data frame
all(targets$Sample == rownames(MDS.data))
# Yes, so we can annotated the coordinates using the 'targets' object
MDS.data = cbind(
  MDS.data, Animal=targets$Animal,
  Infection=targets$Infection,
  Time=targets$Time)
# This way, we preserve the factors defined in the 'targets' object

# Examine the annotated MDS coordinates
head(MDS.data)
# The time factor is still well ordered
MDS.data$Time

# Draw the MDS plot using ggplot2
pdf(file = 'MDS.general.pdf', width = 12, height = 10)
ggplot(MDS.data) +
  geom_point(aes(
    x=X1,
    y=X2,
    colour=Infection,
    shape=Time),
    size=5) +
  ggtitle('2 hours post-infection') +
  xlab('Dimension 1') +
  ylab('Dimension 2') +
  geom_text(aes(label=Animal, x=X1, y=X2, vjust=-1), size=4) +
  theme(
    title = element_text(size = rel(2)),
    axis.text = element_text(size = rel(2))
    )
dev.off()

# Remove the temporary variable
rm(MDS, MDS.data)


# MDS using only 2H samples -----------------------------------------------


# Extract the subset of samples
MDS.data = plotMDS(
  filt_rm_rRNA_norm[,grep('2H', colnames(filt_rm_rRNA_norm))],
  top=nrow(filt_rm_rRNA_norm))

# Extract the coordinates for each sample
MDS.data = data.frame(MDS.data$cmdscale.out)

# Use the sample name to annotate the coordinate with sample information
MDS.data = cbind(MDS.data, strsplit2(rownames(MDS.data), '_'))

# Rename the columns of the data frame for convenience
colnames(MDS.data) = c('x', 'y', 'Animal', 'Infection', 'Time')

# Force the order of Infection levels to have TB in second position
MDS.data$Infection = factor(MDS.data$Infection, levels = c('CN','TB','MB'))

# Draw the MDS plot using ggplot2
pdf(file = 'MDS.2H.pdf', width = 12, height = 10)
ggplot(MDS.data) +
  geom_point(aes(
    x=x,
    y=y,
    colour=Infection,
    shape=Infection),
    size=5) +
  ggtitle('2 hours post-infection') +
  xlab('Dimension 1') +
  ylab('Dimension 2') +
  geom_text(aes(label=Animal, x=x, y=y, vjust=-1), size=4) +
  theme(
    title = element_text(size = rel(2)),
    axis.text = element_text(size = rel(2))
  )
dev.off()


# MDS using only 6H samples -----------------------------------------------


# Extract the subset of samples
MDS.data = plotMDS(
  filt_rm_rRNA_norm[,grep('6H', colnames(filt_rm_rRNA_norm))],
  top=nrow(filt_rm_rRNA_norm))

# Extract the coordinates for each sample
MDS.data = data.frame(MDS.data$cmdscale.out)

# Use the sample name to annotate the coordinate with sample information
MDS.data = cbind(MDS.data, strsplit2(rownames(MDS.data), '_'))

# Rename the columns of the data frame for convenience
colnames(MDS.data) = c('x', 'y', 'Animal', 'Infection', 'Time')

# Force the order of Infection levels to have TB in second position
MDS.data$Infection = factor(MDS.data$Infection, levels = c('CN','TB','MB'))

# Draw the MDS plot using ggplot2
pdf(file = 'MDS.6H.pdf', width = 12, height = 10)
ggplot(MDS.data) +
  geom_point(aes(
    x=x,
    y=y,
    colour=Infection,
    shape=Infection),
    size=5) +
  ggtitle('6 hours post-infection') +
  xlab('Dimension 1') +
  ylab('Dimension 2') +
  geom_text(aes(label=Animal, x=x, y=y, vjust=-1), size=4) +
  theme(
    title = element_text(size = rel(2)),
    axis.text = element_text(size = rel(2))
  )
dev.off()


# MDS using only 24H samples -----------------------------------------------


# Extract the subset of samples
MDS.data = plotMDS(
  filt_rm_rRNA_norm[,grep('24H', colnames(filt_rm_rRNA_norm))],
  top=nrow(filt_rm_rRNA_norm))

# Extract the coordinates for each sample
MDS.data = data.frame(MDS.data$cmdscale.out)

# Use the sample name to annotate the coordinate with sample information
MDS.data = cbind(MDS.data, strsplit2(rownames(MDS.data), '_'))

# Rename the columns of the data frame for convenience
colnames(MDS.data) = c('x', 'y', 'Animal', 'Infection', 'Time')

# Force the order of Infection levels to have TB in second position
MDS.data$Infection = factor(MDS.data$Infection, levels = c('CN','TB','MB'))

# Draw the MDS plot using ggplot2
pdf(file = 'MDS.24H.pdf', width = 12, height = 10)
ggplot(MDS.data) +
  geom_point(aes(
    x=x,
    y=y,
    colour=Infection,
    shape=Infection),
    size=5) +
  ggtitle('6 hours post-infection') +
  xlab('Dimension 1') +
  ylab('Dimension 2') +
  geom_text(aes(label=Animal, x=x, y=y, vjust=-1), size=4) +
  theme(
    title = element_text(size = rel(2)),
    axis.text = element_text(size = rel(2))
  )
dev.off()


# MDS using only 48H samples -----------------------------------------------


# Extract the subset of samples
MDS.data = plotMDS(
  filt_rm_rRNA_norm[,grep('48H', colnames(filt_rm_rRNA_norm))],
  top=nrow(filt_rm_rRNA_norm))

# Extract the coordinates for each sample
MDS.data = data.frame(MDS.data$cmdscale.out)

# Use the sample name to annotate the coordinate with sample information
MDS.data = cbind(MDS.data, strsplit2(rownames(MDS.data), '_'))

# Rename the columns of the data frame for convenience
colnames(MDS.data) = c('x', 'y', 'Animal', 'Infection', 'Time')

# Force the order of Infection levels to have TB in second position
MDS.data$Infection = factor(MDS.data$Infection, levels = c('CN','TB','MB'))

# Draw the MDS plot using ggplot2
pdf(file = 'MDS.48H.pdf', width = 12, height = 10)
ggplot(MDS.data) +
  geom_point(aes(
    x=x,
    y=y,
    colour=Infection,
    shape=Infection),
    size=5) +
  ggtitle('6 hours post-infection') +
  xlab('Dimension 1') +
  ylab('Dimension 2') +
  geom_text(aes(label=Animal, x=x, y=y, vjust=-1), size=4) +
  theme(
    title = element_text(size = rel(2)),
    axis.text = element_text(size = rel(2))
  )
dev.off()