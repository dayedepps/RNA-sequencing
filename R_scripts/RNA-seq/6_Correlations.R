
# Assumes script 2_MergeDEtables.R was run previously


# Libraries and dependencies ----------------------------------------------


library(ggplot2)


# General set up ----------------------------------------------------------


# Folder where the results of the analysis will be saved
rootDir = 'C:/Users/krue/Documents/Kevin-Logs/2013-08-27_AlvMac_RNAseq/04_complete_analysis'

setwd(rootDir)

# Create a folder to store the correlation plots
correlation.dir = 'logFC_correlation'
dir.create(file.path(rootDir, correlation.dir))


# Plot and estimate the correlation of log fold change values -------------


# Load the DE table containing log2 fold-change values
load('all.contrasts.ensembl71.RData')

# Prepare the maximum range of logFC values for the axes
# (to scale all correlation plots identically)
logFC.range = range(c(
  all.contrasts.annotated$logFC.MB.CN.2H,
  all.contrasts.annotated$logFC.MB.CN.6H,
  all.contrasts.annotated$logFC.MB.CN.24H,
  all.contrasts.annotated$logFC.MB.CN.48H,
  all.contrasts.annotated$logFC.TB.CN.2H,
  all.contrasts.annotated$logFC.TB.CN.6H,
  all.contrasts.annotated$logFC.TB.CN.24H,
  all.contrasts.annotated$logFC.TB.CN.48H,
  all.contrasts.annotated$logFC.MB.TB.2H,
  all.contrasts.annotated$logFC.MB.TB.6H,
  all.contrasts.annotated$logFC.MB.TB.24H,
  all.contrasts.annotated$logFC.MB.TB.48H)) +
  c(-1, 1)

# Prepare a mapping between flag and text label
sig.flag = c('None', 'M. tuberculosis', 'M.bovis', 'Both')



# 2H
# Extract the log2FC of the appropriate conditions
cor.data = data.frame(
  TB=all.contrasts.annotated$logFC.TB.CN.2H,
  MB=all.contrasts.annotated$logFC.MB.CN.2H,
  DE=factor(
    sig.flag[
      1+(all.contrasts.annotated$FDR.TB.CN.2H < 0.05) +
        2*(all.contrasts.annotated$FDR.MB.CN.2H < 0.05)
      ],
    levels=sig.flag)
)
# Generate the figure
pdf(file = file.path(correlation.dir, 'correlation_2h.pdf'), width = 8, height = 6)
ggplot(cor.data) +
  geom_point(aes(x=TB, y=MB, colour=DE, shape=DE)) +
  geom_abline(intercept = 0, slope = 1) +
  scale_x_continuous(limits=logFC.range) +
  scale_y_continuous(limits=logFC.range) +
  scale_colour_discrete('Differentially\nexpressed', drop = FALSE) +
  scale_shape_discrete('Differentially\nexpressed', drop = FALSE) +
  ggtitle('2H') +
  xlab('M. tuberculosis - control') +
  ylab('M. bovis - control')
dev.off()

# Intercept, slope, and correlation coefficient of linear model
summary(lm(MB ~ TB, cor.data))
# Remove the temporary object
rm(cor.data)

# 6H
# Extract the log2FC of the appropriate conditions
cor.data = data.frame(
  TB=all.contrasts.annotated$logFC.TB.CN.6H,
  MB=all.contrasts.annotated$logFC.MB.CN.6H,
  DE=factor(
    sig.flag[
      1+(all.contrasts.annotated$FDR.TB.CN.6H < 0.05) +
        2*(all.contrasts.annotated$FDR.MB.CN.6H < 0.05)
      ],
    levels=sig.flag)
)
# Generate the figure
pdf(file = file.path(correlation.dir, 'correlation_6h.pdf'), width = 8, height = 6)
ggplot(cor.data) +
  geom_point(aes(x=TB, y=MB, colour=DE, shape=DE)) +
  geom_abline(intercept = 0, slope = 1) +
  scale_x_continuous(limits=logFC.range) +
  scale_y_continuous(limits=logFC.range) +
  scale_colour_discrete('Differentially\nexpressed', drop = FALSE) +
  scale_shape_discrete('Differentially\nexpressed', drop = FALSE) +
  ggtitle('6H') +
  xlab('M. tuberculosis - control') +
  ylab('M. bovis - control')
dev.off()

# Intercept, slope, and correlation coefficient of linear model
summary(lm(MB ~ TB, cor.data))
# Remove the temporary object
rm(cor.data)

# 24H
# Extract the log2FC of the appropriate conditions
cor.data = data.frame(
  TB=all.contrasts.annotated$logFC.TB.CN.24H,
  MB=all.contrasts.annotated$logFC.MB.CN.24H,
  DE=factor(
    sig.flag[
      1+(all.contrasts.annotated$FDR.TB.CN.24H < 0.05) +
        2*(all.contrasts.annotated$FDR.MB.CN.24H < 0.05)
      ],
    levels=sig.flag)
)
# Generate the figure
pdf(file = file.path(correlation.dir, 'correlation_24h.pdf'), width = 8, height = 6)
ggplot(cor.data) +
  geom_point(aes(x=TB, y=MB, colour=DE, shape=DE)) +
  geom_abline(intercept = 0, slope = 1) +
  scale_x_continuous(limits=logFC.range) +
  scale_y_continuous(limits=logFC.range) +
  scale_colour_discrete('Differentially\nexpressed', drop = FALSE) +
  scale_shape_discrete('Differentially\nexpressed', drop = FALSE) +
  ggtitle('24H') +
  xlab('M. tuberculosis - control') +
  ylab('M. bovis - control')
dev.off()

# Intercept, slope, and correlation coefficient of linear model
summary(lm(MB ~ TB, cor.data))
# Remove the temporary object
rm(cor.data)

# 48H
# Extract the log2FC of the appropriate conditions
cor.data = data.frame(
  TB=all.contrasts.annotated$logFC.TB.CN.48H,
  MB=all.contrasts.annotated$logFC.MB.CN.48H,
  DE=factor(
    sig.flag[
      1+(all.contrasts.annotated$FDR.TB.CN.48H < 0.05) +
        2*(all.contrasts.annotated$FDR.MB.CN.48H < 0.05)
      ],
    levels=sig.flag)
)
# Generate the figure
pdf(file = file.path(correlation.dir, 'correlation_48h.pdf'), width = 8, height = 6)
ggplot(cor.data) +
  geom_point(aes(x=TB, y=MB, colour=DE, shape=DE)) +
  geom_abline(intercept = 0, slope = 1) +
  scale_x_continuous(limits=logFC.range) +
  scale_y_continuous(limits=logFC.range) +
  scale_colour_discrete('Differentially\nexpressed', drop = FALSE) +
  scale_shape_discrete('Differentially\nexpressed', drop = FALSE) +
  ggtitle('48H') +
  xlab('M. tuberculosis - control') +
  ylab('M. bovis - control')
dev.off()

# Intercept, slope, and correlation coefficient of linear model
summary(lm(MB ~ TB, cor.data))
# Remove the temporary object
rm(cor.data)

# Cleanup the temporary objects
rm(all.contrasts.annotated, correlation.dir, logFC.range, sig.flag)
