#############################################################################################################
##Author :          John S. House - BRC;NCSU
##Last Modified :   8.18.18 - John House
##Purpose:          Use Cleaned hash and count file to create normalized counts and max dose contrasts
##NOTE:             ***THIS SCRIPT TAKES THE LONGEST TIME AND THE MOST COMPUTING RESOURCES***
##NOTE:             18 workers at 3.0ghz will take ~ 2-6 minutes on a 3000 feature x 200 sample array for each
##NOTE:               DESeq2 call including each contrast for this 2800 gene dataset
##Info:             MRO version 3.5
##Input:            Hash file and count matrix from TOS_S1....R
##Output:           Normalized Count Matrix (Output Files/Cleaned_data/normcounts.csv) for entire experiment
##Output:           Full model DDS (/Deseq2 dds Files) and each contrast dds (/MaxDoseContrasts)
##NOTE:             This script addresses changes made to DESeq2() that necessitated a re-writing of how
##                  to calculate DEG expresssion, and also a mistake in the initial code that defined 
##                  the grouping variable to include all doses for a chemical vs. the Vehicle instead of 
##                  only the "max" dose of 10.
#############################################################################################################
package.list <- c("DESeq2", "BiocParallel", "dplyr", "ggplot2")
if (length(setdiff(package.list, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(package.list, rownames(installed.packages())))
}
lapply(package.list, require, character.only = TRUE)
setwd("C:/Users/jshou/Google Drive/TempOSeq_Pipeline/ht_cba_manuscript")

load("Output Files/Cleaned_data/qc_counts.RData")
load("Output Files/Cleaned_data/qc_hash.RData")

qc_hash$group <- paste0(qc_hash$Treatment, "_", round(qc_hash$numeric_dose,0))
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

### This is needed for the following function ###
results_dds <- results(dds)

DESeq2_sContrast <- function(deseq_ds = dds, trt, ref, trt_label=trt, ref_label=ref, results_ds = results_dds) {
  testshrunk<-lfcShrink(deseq_ds, contrast = c("group",trt,ref), res = results_ds) ## we want shrunken l2fc ###
  test <- results(deseq_ds, contrast = c("group",trt,ref),parallel = T)  ### everything else is from un-shrunken ###
  test@listData$log2FoldChange <- testshrunk@listData$log2FoldChange  ### replace L2fc estimates with shrunken ones ###
  testdf<-as.data.frame(test);  total<- sum(testdf$padj<.1,na.rm=T)
  testdf$trt <- paste0(trt,"_vs_",ref)
  testdf$Treatment <- paste0(trt_label, " vs. ", ref_label)
  testdf$gene<-rownames(testdf); testdf$id<-paste0(testdf$gene,"/",testdf$Treatment)
  rownames(testdf)<-NULL
  write.csv(testdf,file = paste0("Output Files/MaxDoseContrasts/", trt,"_",ref,".csv"))
  p <- plotMA(test,alpha=0.10, main=paste0(trt_label, " vs. ", ref_label," DEGs = ",total), cex=1,ylim=c(-4,4),cex.main=.9)
  print(p)
  return(testdf)
}
## This function will:
# 1) Create a MA Plot showing the Differentially expressed genes for the comparions
# 2) Use the lfcShrink() function to get the l2fc values
# 3) Use the results() function to get pvalues, adj.pvalues etc..
# 4) Combined into a dataset for all genes tested for that comparison and write a .csv
# 5) return that dataset to an object in memory (i.e. - Dof_10_DEG)
Dof_10_DEG <- DESeq2_sContrast(deseq_ds = dds,
                               trt = "Dofetilide_10",
                               ref = "VEHICLE_0",
                               trt_label = "Max.Dose Dofetilide",
                               ref_label = "Control",
                               results_ds = results_dds)

### feel free to put a wrapper around it for printing PDF's of the graphs ###
pdf(file = "test.pdf", height = 10, width = 10)
par(mfrow = c(2,2)) #set graphs to be a 2 by 2 grid
myds1 <- DESeq2_sContrast( trt = "Dofetilide_10",
                  ref = "VEHICLE_0") 
myds2 <- DESeq2_sContrast( trt = "Isoproterenol_10",
                  ref = "VEHICLE_0")
myds3 <- DESeq2_sContrast( trt = "Nifedipine_10",
                  ref = "VEHICLE_0")
myds4 <- DESeq2_sContrast( trt = "Propranolol_10",
                  ref = "VEHICLE_0")
dev.off()

combined_maxdose <- rbind(myds1, myds2, myds3, myds4)
combined_maxdose$Treatment <- factor(combined_maxdose$Treatment)
save(combined_maxdose, file = "Output Files/MaxDoseContrasts/combined_maxdose.RData")
write.csv(combined_maxdose, file = "Output Files/MaxDoseContrasts/combined_maxdose.csv")


