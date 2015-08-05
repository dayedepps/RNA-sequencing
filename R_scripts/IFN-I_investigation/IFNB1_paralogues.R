
# Assumes script 3_DElists.R was run previously


# General set up ----------------------------------------------------------


# Folder where the results of the analysis will be saved
rootDir = 'C:/Users/krue/Documents/Kevin-Logs/2013-08-27_AlvMac_RNAseq/04_complete_analysis'

setwd(rootDir)

# Create a folder to store raw counts
IFN1.folder = 'IFN-I'
dir.create(file.path(rootDir, IFN1.folder))


# Extract raw counts of IFNB1 paralogues ----------------------------------

# Load the raw counts
load('raw.RData')

# Import the list of IFNB1 paralogues exported from Ensembl 71
IFNB1_paralogs = scan(
  'C:/Users/krue/Documents/Kevin-Logs/2013-08-27_AlvMac_RNAseq/05_IFNa-b/Ensembl71_paralogues/Paralog-71-ENSBTAG00000047142-Ensembl71.txt',
  what = character())

# Fetch the raw counts for each of those Ensembl identifiers
raw.counts = t(raw$counts[IFNB1_paralogs,])

# Export those counts to a csv file
write.table(
  x = raw.counts,
  file = file.path(IFN1.folder, 'raw_counts.txt'),
  quote = FALSE,
  sep = '\t')
