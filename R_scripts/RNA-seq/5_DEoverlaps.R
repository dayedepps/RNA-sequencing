
# Assumes script 3_DElists.R was run previously


# Libraries and dependencies ----------------------------------------------


library(VennDiagram)


# General set up ----------------------------------------------------------


# Folder where the results of the analysis will be saved
rootDir = 'C:/Users/krue/Documents/Kevin-Logs/2013-08-27_AlvMac_RNAseq/04_complete_analysis'

setwd(rootDir)


# Draw 3-way Venn diagrams to display overlaps ----------------------------


# Load the lists of DE genes in each contrasts
load('DElists.RData') 

# Create a folder to save the files containing the Venn diagrams
venn.path = file.path(rootDir, 'Venns')
dir.create(venn.path)

# 2 hours
venn = venn.diagram(
  cex=2,
  fill = c('green', 'yellow', 'darkorchid1'),
  x=DElists[c('MB.CN.2','TB.CN.2','MB.TB.2')],
  filename=NULL)
# Save figure to file
pdf(file = file.path(venn.path, 'Venn_3way_2h.pdf'), width = 10, height = 10)
grid.newpage()
grid.draw(venn)
dev.off()
# Remove the venn object
rm(venn)

# 6 hours
venn = venn.diagram(
  cex=2,
  fill = c('green', 'yellow', 'darkorchid1'),
  x=DElists[c('MB.CN.6','TB.CN.6','MB.TB.6')],
  filename=NULL)
# Save figure to file
pdf(file = file.path(venn.path, 'Venn_3way_6h.pdf'), width = 10, height = 10)
grid.newpage()
grid.draw(venn)
dev.off()
# Remove the venn object
rm(venn)

# 24 hours
venn = venn.diagram(
  cex=2,
  fill = c('green', 'yellow', 'darkorchid1'),
  x=DElists[c('MB.CN.24','TB.CN.24','MB.TB.24')],
  filename=NULL)
# Save figure to file
pdf(file = file.path(venn.path, 'Venn_3way_24h.pdf'), width = 10, height = 10)
grid.newpage()
grid.draw(venn)
dev.off()
# Remove the venn object
rm(venn)

# 48 hours
venn = venn.diagram(
  cex=2,
  fill = c('green', 'yellow', 'darkorchid1'),
  x=DElists[c('MB.CN.48','TB.CN.48','MB.TB.48')],
  filename=NULL)
# Save figure to file
pdf(file = file.path(venn.path, 'Venn_3way_48h.pdf'), width = 10, height = 10)
grid.newpage()
grid.draw(venn)
dev.off()
# Remove the venn object
rm(venn)


# Draw 2-way Venn diagrams to display overlaps ----------------------------


# 2 hours
venn = venn.diagram(
  cex=2,
  fill = c('green', 'yellow'),
  x=DElists[c('MB.CN.2','TB.CN.2')],
  filename=NULL)
# Save figure to file
pdf(file = file.path(venn.path, 'Venn_2way_2h.pdf'), width = 10, height = 10)
grid.newpage()
grid.draw(venn)
dev.off()
# Remove the venn object
rm(venn)

# 6 hours
venn = venn.diagram(
  cex=2,
  fill = c('green', 'yellow'),
  x=DElists[c('MB.CN.6','TB.CN.6')],
  filename=NULL)
# Save figure to file
pdf(file = file.path(venn.path, 'Venn_2way_6h.pdf'), width = 10, height = 10)
grid.newpage()
grid.draw(venn)
dev.off()
# Remove the venn object
rm(venn)

# 24 hours
venn = venn.diagram(
  cex=2,
  fill = c('green', 'yellow'),
  x=DElists[c('MB.CN.24','TB.CN.24')],
  filename=NULL)
# Save figure to file
pdf(file = file.path(venn.path, 'Venn_2way_24h.pdf'), width = 10, height = 10)
grid.newpage()
grid.draw(venn)
dev.off()
# Remove the venn object
rm(venn)

# 48 hours
venn = venn.diagram(
  cex=2,
  fill = c('green', 'yellow'),
  x=DElists[c('MB.CN.48','TB.CN.48')],
  filename=NULL)
# Save figure to file
pdf(file = file.path(venn.path, 'Venn_2way_48h.pdf'), width = 10, height = 10)
grid.newpage()
grid.draw(venn)
dev.off()
# Remove the venn object
rm(venn)


# Cleanup temporary variables ---------------------------------------------


# Remove the path to venn folder and DE lists
rm(venn.path, DElists)
