
# Assumes script 3_DElists.R was run previously


# Libraries and dependencies ----------------------------------------------


library(biomaRt)
library(sigora)


# General set up ----------------------------------------------------------


# Folder where the results of the analysis will be saved
rootDir = 'C:/Users/krue/Documents/Kevin-Logs/2013-08-27_AlvMac_RNAseq/04_complete_analysis'

setwd(rootDir)


# Functions to automate the analyses --------------------------------------

# This function converts bovine Ensembl gene identifiers to human ones
# Bovine genes mapping to a unique human genes are directly converted
# Bovine genes mapping to no human gene are ignored
# For bovine genes mapping to multiple human genes,
# gene names are compared
# if the bovine gene name matches the human gene name for one of the multiple
# hits, the first human gene identifier is used
# otherwise, the bovine gene is ignored
bovine2human_ensembl = function(bovine_ensembl){
  # Connect to bovine Ensembl mart release 71
  mart_bovine = useMart(
    host = 'apr2013.archive.ensembl.org',
    biomart='ENSEMBL_MART_ENSEMBL', dataset='btaurus_gene_ensembl')
  # Connect to human Ensembl mart release 71
  mart_human = useMart(
    host = 'apr2013.archive.ensembl.org',
    biomart = 'ENSEMBL_MART_ENSEMBL', dataset = 'hsapiens_gene_ensembl')
  # get the human_ortholog_ensembl and bovine_gene_symbol
  BM1 = getBM(
    attributes=c(
      'ensembl_gene_id', 'external_gene_id', 'hsapiens_homolog_ensembl_gene'
      ),
    filters='ensembl_gene_id',
    values=bovine_ensembl,
    mart=mart_bovine)
  names(BM1) = c('bovine_ensembl','bovine_symbol','human_ensembl')
  # Get the human_gene_symbol
  BM2 = getBM(
    attributes=c('ensembl_gene_id', 'external_gene_id'),
    filters='ensembl_gene_id',
    values=BM1$human_ensembl,
    mart=mart_human
    )
  names(BM2) = c('human_ensembl', 'human_symbol')
  # Merge the two tables, keeping only the rows mapped to a human_ensembl
  # (bovine_ensembl without human_ensembl cannot be taken any further)
  BM = merge(x=BM1, y=BM2, by.x='human_ensembl')
  # Summarise number of matches for each bovine ensembl
  BM.matches = aggregate(
    formula=human_ensembl~bovine_ensembl,
    data=BM,
    FUN=length
    )
  # Take the human_ensembl where the bovine_ensembl had a unique match
  human_ensembl = BM$human_ensembl[
    which(
      BM$bovine_ensembl %in% BM.matches$bovine_ensembl[BM.matches$human_ensembl == 1]
      )
    ]
  # Remove the unique matches found above
  BM.2more = BM[!BM$bovine_ensembl %in% BM.matches$bovine_ensembl[BM.matches$human_ensembl == 1],]
  # Prepares a variable to store human ensembl mapped to bovine through their
  # gene symbol after being found multimapped on their ensembl id
  human_repicked = c()
  if(nrow(BM.2more) != 0){
    # Identify the rows where the symbol match
    BM.2more$symbol.match = BM.2more$bovine_symbol == BM.2more$human_symbol
    BM.symbol.matches = aggregate(
      formula=symbol.match~bovine_ensembl,
      data=BM.2more,
      FUN=sum
      )
    # Ignore the ones with 0 symbol matches (even though they have multiple
    # ensembl matches)
    # (These have either no bovine symbol,
    # or symbol different from all human symbols
    # For the ones with a unique symbol match, take the corresponding
    # human_ensembl
    # (These are the best match)
    for(bov_ens in BM.symbol.matches$bovine_ensembl[BM.symbol.matches$symbol.match == 1]){
      human_repicked = c(
        human_repicked,
        BM.2more$human_ensembl[
          BM.2more$bovine_ensembl == bov_ens &
            BM.2more$symbol.match == TRUE
          ]
        )
      rm(bov_ens)
    }
    # For the ones with multiple symbols matches (case not seen),
    # take the first human_ensembl
    # (Same human symbol for all, therefore likely to be annotated to the same
    # pathway if annnotated at all)
    # Also, only one human_ensembl to replace one bovine_ensembl, as more
    # than one could skew the enrichment.
    for(bov_ens in BM.symbol.matches$bovine_ensembl[BM.symbol.matches$symbol.match >= 2]){
      human_repicked = c(
        human_repicked,
        BM.2more$human_ensembl[BM.2more$bovine_ensembl == bov_ens][1]
        )
      rm(bov_ens)
    }
  }
  # Merge all the human_ensembl found at each step
  human_ensembl.final = c(human_ensembl, human_repicked)
  # return the result
  return(human_ensembl.final)
}


# This function takes a vector of bovine Ensembl gene identifiers and
# saves the results of a SIGORA analysis after converting them to human
# gene identifiers using the above 'bovine2human_ensembl' function
SIGORA.bovine = function(bovine_ensembl, folder.out, basename.out){
  # convert the bovine ensembl into human ensembl
  human_ensembl = bovine2human_ensembl(bovine_ensembl=bovine_ensembl)
  # convert query list into SIGORA identifiers
  myquerylist <<- ens_converter(x=human_ensembl)
  # Writes the screen output to report file
  sink(
    file=paste(folder.out, '/', basename.out, '_SIGORA.sink.txt', sep=''),
    split=T,
    append=F
    )
  # perform signature over-representation analysis using KEGG GPS
  sigs(samplename=myquerylist, archive='k', markers=1, level=2)
  # print the list of pathways with p-values
  print(summary_results)
  # stop piping into the file
  sink()
  # Painful to read, but contains all the useful information
  export_results(
    filename=paste(folder.out, '/', basename.out, '_SIGORA.results.txt', sep='')
    )
  # Cleanup
  Sys.unsetenv(myquerylist)
}


# Running SIGORA analyses using ensembl 71 annotations --------------------


# Load the lists of DE genes in each contrasts
load('DElists.RData')

# Create a folder to save the files containing differential expression results
sigora71.path = file.path(rootDir, 'SIGORA_Ensembl71')
dir.create(sigora71.path)

# 2 hours post-infection
# MB-CN genes (90)
MB.CN.2 = DElists$MB.CN.2
SIGORA.bovine(bovine_ensembl=MB.CN.2, folder.out=sigora71.path, basename.out='MB.CN.2')
rm(MB.CN.2)
# TB-CN genes (15)
TB.CN.2 = DElists$TB.CN.2
SIGORA.bovine(bovine_ensembl=TB.CN.2, folder.out=sigora71.path, basename.out='TB.CN.2')
rm(TB.CN.2)
# 75 MB-specific genes at 2h
MB.spe.2 = DElists$MB.CN.2[!DElists$MB.CN.2 %in% DElists$TB.CN.2]
SIGORA.bovine(bovine_ensembl=MB.spe.2, folder.out=sigora71.path, basename.out='MB.spe.2')
rm(MB.spe.2)
# No TB-specific DE gene at 2h
# No MB-TB DE genes at 2h
# Only 15 DE genes common to MB and TB at 2h
common.2 = DElists$MB.CN.2 [DElists$MB.CN.2 %in% DElists$TB.CN.2]
SIGORA.bovine(bovine_ensembl=common.2, folder.out=sigora71.path, basename.out='common.2')
# Quick look at the genes and gene description because there are only 15
load('topTags/all.contrasts.annotated.RData')
all.contrasts.annotated[all.contrasts.annotated$ensembl_gene_id %in% common.2,c('external_gene_id','description')]
rm(all.contrasts.annotated.RData)
rm(common.2)

# 6 hours post-infection
# MB-CN genes (1386)
MB.CN.6 = DElists$MB.CN.6
SIGORA.bovine(bovine_ensembl=MB.CN.6, folder.out=sigora71.path, basename.out='MB.CN.6')
rm(MB.CN.6)
# TB-CN genes (1969)
TB.CN.6 = DElists$TB.CN.6
SIGORA.bovine(bovine_ensembl=TB.CN.6, folder.out=sigora71.path, basename.out='TB.CN.6')
rm(TB.CN.6)
# 381 MB-specific DE genes at 6h
MB.spe.6 = DElists$MB.CN.6 [!DElists$MB.CN.6 %in% DElists$TB.CN.6]
SIGORA.bovine(bovine_ensembl=MB.spe.6, folder.out=sigora71.path, basename.out='MB.spe.6')
rm(MB.spe.6)
# 964 TB-specific genes at 6h
TB.spe.6 = DElists$TB.CN.6 [!DElists$TB.CN.6 %in% DElists$MB.CN.6]
SIGORA.bovine(bovine_ensembl=TB.spe.6, folder.out=sigora71.path, basename.out='TB.spe.6')
rm(TB.spe.6)
# No MB-TB DE genes at 6h
# 1005 genes common to MB and TB at 6h
common.6 = DElists$MB.CN.6 [DElists$MB.CN.6 %in% DElists$TB.CN.6]
SIGORA.bovine(bovine_ensembl=common.6, folder.out=sigora71.path, basename.out='common.6')
rm(common.6)

# 24 hours post-infection
# MB-CN genes (5726)
MB.CN.24 = DElists$MB.CN.24
SIGORA.bovine(bovine_ensembl=MB.CN.24, folder.out=sigora71.path, basename.out='MB.CN.24')
rm(MB.CN.24)
# TB-CN genes (4224)
TB.CN.24 = DElists$TB.CN.24
SIGORA.bovine(bovine_ensembl=TB.CN.24, folder.out=sigora71.path, basename.out='TB.CN.24')
rm(TB.CN.24)
# 1721 MB-specific DE genes at 24h
MB.spe.24 = DElists$MB.CN.24 [!DElists$MB.CN.24 %in% DElists$TB.CN.24]
SIGORA.bovine(bovine_ensembl=MB.spe.24, folder.out=sigora71.path, basename.out='MB.spe.24')
rm(MB.spe.24)
# 219 TB-specific genes at 24h
TB.spe.24 = DElists$TB.CN.24 [!DElists$TB.CN.24 %in% DElists$MB.CN.24]
SIGORA.bovine(bovine_ensembl=TB.spe.24, folder.out=sigora71.path, basename.out='TB.spe.24')
rm(TB.spe.24)
# No MB-TB DE genes at 24h
# 4005 genes common to MB and TB at 24h
common.24 = DElists$MB.CN.24 [DElists$MB.CN.24 %in% DElists$TB.CN.24]
SIGORA.bovine(bovine_ensembl=common.24, folder.out=sigora71.path, basename.out='common.24')
rm(common.24)

# 48 hours post-infection
# MB-CN genes (7514)
MB.CN.48 = DElists$MB.CN.48
SIGORA.bovine(bovine_ensembl=MB.CN.48, folder.out=sigora71.path, basename.out='MB.CN.48')
rm(MB.CN.48)
# TB-CN genes (4484)
TB.CN.48 = DElists$TB.CN.48
SIGORA.bovine(bovine_ensembl=TB.CN.48, folder.out=sigora71.path, basename.out='TB.CN.48')
rm(TB.CN.48)
# 3441 MB-specific DE genes at 48h
MB.spe.48 = DElists$MB.CN.48 [!DElists$MB.CN.48 %in% DElists$TB.CN.48]
SIGORA.bovine(bovine_ensembl=MB.spe.48, folder.out=sigora71.path, basename.out='MB.spe.48')
rm(MB.spe.48)
# 411 TB-specific genes at 48h
TB.spe.48 = DElists$TB.CN.48 [!DElists$TB.CN.48 %in% DElists$MB.CN.48]
SIGORA.bovine(bovine_ensembl=TB.spe.48, folder.out=sigora71.path, basename.out='TB.spe.48')
rm(TB.spe.48)
# 3438 MB-TB DE genes at 48h
diff.48 = DElists$MB.TB.48
SIGORA.bovine(bovine_ensembl=diff.48, folder.out=sigora71.path, basename.out='MB.TB.diff.48')
rm(diff.48)
# 4073 genes common to MB and TB at 48h
common.48 = DElists$MB.CN.48 [DElists$MB.CN.48 %in% DElists$TB.CN.48]
SIGORA.bovine(bovine_ensembl=common.48, folder.out=sigora71.path, basename.out='common.48')
rm(common.48)
# 2565 genes no different between MB and TB, but different between each and CN
tmp = intersect(DElists$MB.CN.48, DElists$TB.CN.48)
similar.48 = tmp [which(!tmp %in% DElists$MB.TB.48)]
SIGORA.bovine(bovine_ensembl=similar.48, folder.out=sigora71.path, basename.out='similar.48')
rm(tmp, similar.48)
# DE genes different between all three infection groups (CN, TB, MB)
tmp = intersect(DElists$MB.CN.48, DElists$TB.CN.48)
all.different.48 = intersect(tmp, DElists$MB.TB.48)
SIGORA.bovine(bovine_ensembl=all.different.48, folder.out=sigora71.path, basename.out='all.different.48')
rm(tmp, all.different.48)


rm(DElists, myquerylist, sigora71.path, SIGORA.bovine, bovine2human_ensembl)
