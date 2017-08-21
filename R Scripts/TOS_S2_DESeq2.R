#############################################################################################################
##Author :          John S. House - BRC;NCSU
##Last Modified :   8.16.17 - John House 
##Purpose:          Use Cleaned hash and count file to create normalized counts and max dose contrasts
##NOTE:             ***THIS SCRIPT TAKES THE LONGEST TIME AND THE MOST COMPUTING RESOURCES***
##NOTE:             18 workers at 3.0ghz will take ~ 2-6 minutes on a 3000 feature x 200 sample array for each
##NOTE:               DESeq2 call including each contrast
##Info:             MRO version 3.4
##Input:            Hash file and count matrix from TOS_S1....R
##Output:           Normalized Count Matrix (Output Files/Cleaned_data/normcounts.csv) for entire experiment
##Output:           Full model DDS (/Deseq2 dds Files) and each contrast dds (/MaxDoseContrasts)
#############################################################################################################
package.list <- c("DESeq2", "BiocParallel", "dplyr", "ggplot2")
if (length(setdiff(package.list, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(package.list, rownames(installed.packages())))
}
lapply(package.list, require, character.only = TRUE)
setwd("C:/Users/jshouse/Google_Drive/TempOSeq_Pipeline/ht_cba_manuscript")

load("Output Files/Cleaned_data/qc_counts.RData")
load("Output Files/Cleaned_data/qc_hash.RData")

qc_hash$group <- factor(qc_hash$group)

## group was defined as one level up from the lowest (replicate) ##
## you can define another level to test if you like ##
## this experimental design still assumes that there are different media controls on each Index.Set ##
dds = DESeqDataSetFromMatrix(countData = qc_counts,
                             colData = qc_hash,
                             design = ~ group)

system.time(dds <- DESeq(dds, parallel = TRUE))
save(dds, file = "Output Files/DESeq2 dds Files/full_dds.RData")
normcounts <- counts(dds, normalized = TRUE)
write.csv(normcounts, file = "Output Files/Cleaned_data/normcounts.csv")
save(normcounts, file = "Output Files/Cleaned_data/normcounts.RData")

### WRITING OUT ALL INDIVIDUAL CONTRASTS ###
### IN THIS EXPERIMENT THE MAX DOSE WAS 10 mM; THIS WILL NEED TO CHANGE FOR OTHER EXPERIMENTS###
colData(dds)
contrast_hash <- as.data.frame(colData(dds))
treatment <-
  select(dplyr::filter(contrast_hash, DOSE == 10),
         group,
         Treatment,
         DOSE)
treatment$group <- as.character(treatment$group)
treatment$control <- "VEHICLE"

### INITIALIZE EMPTY DATAFRAME ###
combined_maxdose <- data.frame()
ptm <- proc.time()
for (i in 1:nrow(treatment)) {
  print(paste0("Starting contrast ", i, " of ", nrow(treatment)))
  myds <- treatment[i, ]
  test <-
    results(
      dds,
      contrast = c("group", myds$group[1], myds$control[1]),
      parallel = TRUE
    )
  save(test,
       file = paste0("Output Files/MaxDoseContrasts/", myds$Treatment[1], ".RData"))
  testdf <- as.data.frame(test)
  testdf$Treatment <- myds$Treatment[1]
  testdf$gene <-
    rownames(testdf)
  testdf$id <- paste0(testdf$gene, "/", testdf$Treatment)
  rownames(testdf) <- NULL
  write.csv(testdf,
            file = paste0("Output Files/MaxDoseContrasts/", myds$Treatment[1], ".csv"))
  combined_maxdose <- rbind(combined_maxdose, testdf)
}
proc.time() - ptm
combined_maxdose$Treatment <- factor(combined_maxdose$Treatment)
save(combined_maxdose, file = "Output Files/MaxDoseContrasts/combined_maxdose.RData")
write.csv(combined_maxdose, file = "Output Files/MaxDoseContrasts/combined_maxdose.csv")




