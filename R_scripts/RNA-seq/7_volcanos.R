
# Assumes script 2_MergeDEtables.R was run previously


# Libraries and dependencies ----------------------------------------------


library(ggplot2)


# General set up ----------------------------------------------------------


# Folder where the results of the analysis will be saved
rootDir = 'C:/Users/krue/Documents/Kevin-Logs/2013-08-27_AlvMac_RNAseq/04_complete_analysis'

setwd(rootDir)

# Create a folder to store the volcano plots
volcano.dir = 'volcanos'
dir.create(file.path(rootDir, volcano.dir))


# Plot the p-value against the fold-change --------------------------------


# Import the merged DE table produced previously
load('all.contrasts.ensembl71.RData')

# 2H
volcano.data = data.frame(
  log2FC=all.contrasts.annotated$logFC.MB.TB.2H,
  mlog10FDR=-log10(all.contrasts.annotated$PValue.MB.TB.2H),
  sig=factor(
    all.contrasts.annotated$FDR.MB.TB.2H < 0.05, levels = c(FALSE, TRUE)
    )
  )
# Generate the figure
pdf(file = file.path(volcano.dir, 'volcano_2h.pdf'))
ggplot(volcano.data) +
  geom_point(aes(x=log2FC, y=mlog10FDR, colour=sig)) +
  scale_colour_manual(
    drop = FALSE, guide = "none", values = c('black', 'red')) +
  xlab('log2 fold change') +
  ylab('-log10 FDR') +
  ggtitle('2H') +
  scale_x_continuous(limits=c(-4, 4)) +
  scale_y_continuous(limits=c(0, 15))
dev.off()
# Remove the temporary object
rm(volcano.data)

# 6H
volcano.data = data.frame(
  log2FC=all.contrasts.annotated$logFC.MB.TB.6H,
  mlog10FDR=-log10(all.contrasts.annotated$PValue.MB.TB.6H),
  sig=factor(
    all.contrasts.annotated$FDR.MB.TB.6H < 0.05, levels = c(FALSE, TRUE)
  )
)
# Generate the figure
pdf(file = file.path(volcano.dir, 'volcano_6h.pdf'))
ggplot(volcano.data) +
  geom_point(aes(x=log2FC, y=mlog10FDR, colour=sig)) +
  scale_colour_manual(
    drop = FALSE, guide = "none", values = c('black', 'red')) +
  xlab('log2 fold change') +
  ylab('-log10 FDR') +
  ggtitle('6H') +
  scale_x_continuous(limits=c(-4, 4)) +
  scale_y_continuous(limits=c(0, 15))
dev.off()
# Remove the temporary object
rm(volcano.data)

# 24h
volcano.data = data.frame(
  log2FC=all.contrasts.annotated$logFC.MB.TB.24H,
  mlog10FDR=-log10(all.contrasts.annotated$PValue.MB.TB.24H),
  sig=factor(
    all.contrasts.annotated$FDR.MB.TB.24H < 0.05, levels = c(FALSE, TRUE)
  )
)
# Generate the figure
pdf(file = file.path(volcano.dir, 'volcano_24h.pdf'))
ggplot(volcano.data) +
  geom_point(aes(x=log2FC, y=mlog10FDR, colour=sig)) +
  scale_colour_manual(
    drop = FALSE, guide = "none", values = c('black', 'red')) +
  xlab('log2 fold change') +
  ylab('-log10 FDR') +
  ggtitle('24H') +
  scale_x_continuous(limits=c(-4, 4)) +
  scale_y_continuous(limits=c(0, 15))
dev.off()
# Remove the temporary object
rm(volcano.data)

# 48h
volcano.data = data.frame(
  log2FC=all.contrasts.annotated$logFC.MB.TB.48H,
  mlog10FDR=-log10(all.contrasts.annotated$PValue.MB.TB.48H),
  sig=factor(
    all.contrasts.annotated$FDR.MB.TB.48H < 0.05, levels = c(FALSE, TRUE)
  )
)
# Generate the figure
pdf(file = file.path(volcano.dir, 'volcano_48h.pdf'))
ggplot(volcano.data) +
  geom_point(aes(x=log2FC, y=mlog10FDR, colour=sig)) +
  scale_colour_manual(
    drop = FALSE, guide = "none", values = c('black', 'red')) +
  xlab('log2 fold change') +
  ylab('-log10 FDR') +
  ggtitle('48H') +
  scale_x_continuous(limits=c(-4, 4)) +
  scale_y_continuous(limits=c(0, 15))
dev.off()
# Remove the temporary object
rm(volcano.data)

# Cleanup the temporary objects
rm(all.contrasts.annotated, volcano.dir)
