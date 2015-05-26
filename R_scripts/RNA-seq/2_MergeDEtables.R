# Assumes script 1_alv-mac-edgeR-pipeline.R was run previously

# Libraries and dependencies ----------------------------------------------


library(biomaRt)
library(tools)


# General set up ----------------------------------------------------------


# Folder where the results of the analysis will be saved
rootDir = 'C:/Users/krue/Documents/Kevin-Logs/2013-08-27_AlvMac_RNAseq/04_complete_analysis'

setwd(rootDir)


# Import all DE tables ----------------------------------------------------


topTagsfiles = list.files(path='topTags',pattern='*H.RData$')

for(f in topTagsfiles){
  load(file=file.path('topTags', f))
}
rm(f, topTagsfiles)


# Merge all contrasts at all time points in the same dataset --------------


# Progressively merge all tables MB.CN together
all.contrasts = merge(
  x=MB.CN.2$table, y=MB.CN.6$table,
  by='row.names', suffixes=c('.MB.CN.2H','.MB.CN.6H')
  )

colnames(MB.CN.24$table) = paste(colnames(MB.CN.24), 'MB.CN.24H', sep='.')
all.contrasts = merge(
  x=all.contrasts, y=MB.CN.24$table, by.x='Row.names', by.y='row.names'
  )

colnames(MB.CN.48$table) = paste(colnames(MB.CN.48), 'MB.CN.48H', sep='.')
all.contrasts = merge(
  x=all.contrasts, y=MB.CN.48$table, by.x='Row.names', by.y='row.names'
  )

# Progressively add all tables TB.CN
colnames(TB.CN.2$table) = paste(colnames(TB.CN.2$table), 'TB.CN.2H', sep='.')
all.contrasts = merge(
  x=all.contrasts, y=TB.CN.2$table, by.x='Row.names', by.y='row.names'
  )

colnames(TB.CN.6$table) = paste(colnames(TB.CN.6$table), 'TB.CN.6H', sep='.')
all.contrasts = merge(
  x=all.contrasts, y=TB.CN.6$table, by.x='Row.names', by.y='row.names'
  )

colnames(TB.CN.24$table) = paste(colnames(TB.CN.24$table), 'TB.CN.24H', sep='.')
all.contrasts = merge(
  x=all.contrasts, y=TB.CN.24$table, by.x='Row.names', by.y='row.names'
  )

colnames(TB.CN.48$table) = paste(colnames(TB.CN.48$table), 'TB.CN.48H', sep='.')
all.contrasts = merge(
  x=all.contrasts, y=TB.CN.48$table, by.x='Row.names', by.y='row.names'
  )

# Progressively add all tables MB.TB
colnames(MB.TB.2$table) = paste(colnames(MB.TB.2$table), 'MB.TB.2H', sep='.')
all.contrasts = merge(
  x=all.contrasts, y=MB.TB.2$table, by.x='Row.names', by.y='row.names'
  )

colnames(MB.TB.6$table) = paste(colnames(MB.TB.6$table), 'MB.TB.6H', sep='.')
all.contrasts = merge(
  x=all.contrasts, y=MB.TB.6$table, by.x='Row.names', by.y='row.names'
  )

colnames(MB.TB.24$table) = paste(
  colnames(MB.TB.24$table), 'MB.TB.24H', sep='.')
all.contrasts = merge(
  x=all.contrasts, y=MB.TB.24$table, by.x='Row.names', by.y='row.names'
  )

colnames(MB.TB.48$table) = paste(
  colnames(MB.TB.48$table), 'MB.TB.48H', sep='.')
all.contrasts = merge(
  x=all.contrasts, y=MB.TB.48$table, by.x='Row.names', by.y='row.names'
  )

# Examine the resulting table
all.contrasts[1:5,1:5]

# Restore the row names
rownames(all.contrasts) = all.contrasts$Row.names

# Remove the extra column
all.contrasts$Row.names = NULL

# Save the object to file
save(all.contrasts, file=file.path('topTags', 'all.contrasts.RData'))

# Further compress the file to optimise disk space usage
resaveRdaFiles(file.path('topTags', 'all.contrasts.RData'))

# Export the table to a text file
write.table(
  x=cbind(ID=rownames(all.contrasts), all.contrasts),
  file=file.path('topTags', 'all.contrasts.txt'),
  append=F, quote=F, sep='\t', row.names=F, col.names=TRUE
  )

# Remove the individual DE tables
rm(
  list=grep(
    pattern='[[:upper:]]{2}.[[:upper:]]{2}.[[:digit:]]{1,2}',
    x=ls(), value=TRUE
    )
  )



# Annotate the merged table with Ensembl gene identfiers ------------------


# Connect to Ensembl release 71
mart = useMart(
  host = 'apr2013.archive.ensembl.org',
  biomart='ENSEMBL_MART_ENSEMBL',dataset='btaurus_gene_ensembl')

# Get associated gene name and description for all genes 
annotations = getBM(
  attributes=c('ensembl_gene_id','external_gene_id','description'),
  mart=mart)

# Save the Ensembl annotations
# (in case they are removed from the server)
save(annotations, file=file.path('topTags','annotations.Rdata'))

# Add a column containing the Ensembl identifier in the merged table
all.contrasts$ensembl_gene_id = rownames(all.contrasts)

# Merge the annotations with the data
all.contrasts.annotated = merge(
  x=all.contrasts, y=annotations,
  by='ensembl_gene_id', all.x=T, sort=F)

# Save the object to file
save(all.contrasts.annotated, file='all.contrasts.ensembl71.RData')

# Further compress the file to optimise disk space usage
resaveRdaFiles('all.contrasts.ensembl71.RData')

# Export the table to a text file
write.table(
  x=all.contrasts.annotated, file='all.contrasts.ensembl71.txt',
  quote=F, sep='\t', row.names=F)


# Remove the temporary objects
rm(all.contrasts, all.contrasts.annotated, annotations, mart)
