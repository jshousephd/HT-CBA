#############################################################################################################
##Author :          John S. House - BRC;NCSU
##Last Modified :   8/15/2017 - John House
##License:          MIT - Permission to Use and Modify, cite author
##Purpose:          Bring in Hashfile and Count Matrix and Prepare for Analysis. Assess Scaled Raw Media Counts
##Info:             R Version MRO version 3.4.0
##Input:            Hash file and count matrix from TempOseq CArdiomyocyte data prep
##Output:           cleaned qc_count and qc_hash for DESeq2 in /Cleaned_data and
##Output:           low_count_report.csv and PCA of Raw Counts in /RawCountResults
##Output:           histogram of the mean correlations of Raw Counts in /RawCountResults
#############################################################################################################
package.list <-
  c("tidyverse",
    "stringr",
    "scatterD3",
    "plot3D",
    "PerformanceAnalytics")
if (length(setdiff(package.list, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(package.list, rownames(installed.packages())))
}
lapply(package.list, require, character.only = TRUE)

setwd("C:/Users/jshouse/Google_Drive/TempOSeq_Pipeline/ht_cba_manuscript")

### READ IN COUNT MATRIX AND ATTRIBUTE FILE (HASH FILE) ###
raw_count_data <-
  read.csv(
    file = "Input Files/raw_counts.csv",
    header = T,
    stringsAsFactors = F,
    row.names = 1
  ) %>% as.data.frame()
hash <-
  read.csv(
    file = "Input Files/hash.csv",
    stringsAsFactors = FALSE,
    header = T,
    row.names = 13
  )

# Prepare files for DESeq2 ----------------------------------------------------------------------------------------
### remove flagged rows from hash file as well as selecting only columns we need ###
### Hash file can contain other attributes you may want to keep for more complicated modeling in DESeq2 for example ###
hash1 <-
  subset(
    hash,
    is.na(DoNotUse),
    select = c(
      ID,
      Index.Set,
      X96wp.position,
      Treatment,
      Sample.ID,
      DOSE,
      Replicate
    )
  )
### check for missing values in any column ###
try(if (sum(apply(hash1, 2, anyNA)) > 0)
  stop("You have missing values in Hash table, Please fix"))
### replace all non-standard characters in treatment name with '_'   ( could also use  |[:space:]  in clause to replace spaces###
hash1$Treatment <-
  str_replace_all(hash1$Treatment, "[:punct:]", "_")
### Create lowest level of experiment. This is one level up from replicates ####
### This is also the normal defining level for DESeq2 to do contrasts for DEGs analysis ###
hash1$group <- hash1$Treatment
### create numeric dose for use in linear modeling with DESeq2 ###
hash1$numeric_dose <- as.numeric(as.character(hash1$DOSE))
### Raw counts read in added 'X' to first position of each column name because doesn't like starting with 0 ###
### subset count matrix by the rows left in hash1 AND in the same order! ###
finalcounts <- raw_count_data[, rownames(hash1)]
### check for NA's in count matrix
if (sum(is.na(as.matrix(finalcounts))) > 0) {
  stop("Final Count Matrix has NA, please fix")
}
### confirm our colnames of count matrix match our rownames of coldata matrix ###
if (sum(rownames(hash1) != colnames(finalcounts)) > 0) {
  stop("rownames of hash must equal colnames of counts")
}
### REMOVE TRANSCRIPTS THAT HAVE 0 OR 1 as a rowSum across experiment ###
counts_clean1 <-
  finalcounts[rowSums(finalcounts) > 1, ]
nrow(counts_clean1)
summary(counts_clean1)
# mean(as.matrix(counts_clean1))
# mean(as.matrix(finalcounts))
### Assess Well counts ###
### Summarize Each Column Total ###
hash1$col_totals <- colSums(counts_clean1)
hash1$raw_reads <- ifelse(
  hash1$col_totals < 100000,
  "<100k",
  ifelse(
    hash1$col_totals >= 100000 & hash1$col_totals < 200000,
    "100k-200k",
    ifelse(
      hash1$col_totals >= 200000 & hash1$col_totals < 500000,
      "200k-500k",
      ifelse(
        hash1$col_totals >= 500000 &
          hash1$col_totals < 1000000,
        "500k-1m",
        ">1m"
      )
    )
  )
)
table(hash1$raw_reads)

### Experimental Wells to be removed because of Low Total counts  #######
raw_count_cutoff <- 100000
### For our test data, there are no data to be removed due to insufficient counts ###
### Following File will contain the treatments removed if any ###
write.csv(subset(hash1, hash1$col_totals <= raw_count_cutoff),
          file = "Output Files/RawCountResults/low_count_report.csv")
qc_hash <- subset(hash1, hash1$col_totals > raw_count_cutoff)
qc_counts <- counts_clean1[, rownames(qc_hash)]

## Rows (probes/genes) with zero counts or only 1 count across all treatments are removed. ##
## This can be changed by the user to be more or less stringent as needed ###
qc_counts <- qc_counts[rowSums(qc_counts) > 1, ]
nrow(qc_counts)
qc_hash$Treatment <- factor(qc_hash$Treatment)


#  ----------------------------------------------------------------------------------------------------------------
# Replication Assessment for our VEHICLE Controls -----------------------------------------------------------------

counts <- qc_counts
attrs <- qc_hash
attrs$ID2 <- paste0(attrs$group, "_", attrs$Replicate)
head(attrs)
media_counts <-
  counts[, rownames(subset(attrs, Treatment == "VEHICLE"))]
media_counts <- media_counts[rowSums(media_counts) > 0, ]
### Replacing names of count data with the ID2 variable that matches the rownname of the hash file (attrs) ###
mm <- match(names(media_counts), rownames(attrs))
names(media_counts)[!is.na(mm)] <-
  as.character(attrs$ID2[na.omit(mm)])
media_counts <- media_counts[rowSums(media_counts) > 0, ]
### Since these are Raw Counts we want to normalize by library size first ###
media_counts_lib_size <-
  apply(media_counts, 2, function(x)
    x / sum(x))
pcamatrix <- as.matrix(media_counts_lib_size)

### The main source of 'structure' in the data are the high to low counts by gene ###
### We need to get rid of that structure to assess Principal Components of our Vehicle Counts for Replication###

X <- t(scale(t(pcamatrix)))  ## Scales Treatment Columns ##
p = prcomp(X, center = T)
pc <- unclass(p$rotation)[, 1:3]
p_var = (round(100 * (p$sdev ^ 2 / sum(p$sdev ^ 2)), 1))[1:3]

pdf(file = "Output Files/RawCountResults/PCA_MediaControls_RawCounts.pdf",
    width = 10.5,
    height = 8.0)
scatter3D(
  x = pc[, 1],
  y = pc[, 2],
  z = pc[, 3],
  colvar = NULL,
  cex = 2,
  btw = "b2",
  pch = 16,
  col = "darkorchid",
  ticktype = "detailed",
  cex.main = 1.5,
  cex.lab = 1.25,
  cex.axis = .6,
  main = "Principal Components (Controls)",
  xlab = paste0("PC1 - ", p_var[1], "% Variance"),
  ylab = paste0("PC2 - ", p_var[2], "% Variance"),
  zlab = paste0("PC3 - ", p_var[3], "% Variance")
)
text3D(
  pc[, 1] - .05,
  pc[, 2],
  pc[, 3] - .05,
  labels = rep(c(1:24)),
  add = TRUE,
  cex = 1
)
dev.off()

### If 2D Plots are Desired ###
pdf(file = "Output Files/RawCountResults/PCA_Veh_Controls_Raw2D.pdf",
    width = 8,
    height = 6)
par(mar = c(5, 5, 5, 5))
plot(
  x = pc[, 1],
  y = pc[, 2],
  pch = 16,
  col = "darkorchid",
  cex = 2.25,
  main = "Vehicle Control Principal Components",
  cex.main = 2.25,
  xlab = paste0("PC1 - ", p_var[1], "% Variance"),
  cex.lab = 1.75,
  cex.axis = 1.75,
  ylab = paste0("PC2 - ", p_var[2], "% Variance")
)
text(pc[, 1] - .02, pc[, 2] - .01, labels = rep(c(1:24)))
dev.off()

### Assessment of Correlation Structure of the Controls ###
cor_media <- cor(media_counts_lib_size)
meancors <-
  data.frame(Treatment = rownames(cor_media),
             MeanCorr = rowMeans(cor_media))

pdf(file = "Output Files/RawCountResults/Hist_Vehicle_Correlations2D.pdf",
    width = 8,
    height = 6.0)
par(mar = c(5, 5, 5, 5))
hist(
  meancors$MeanCorr,
  main = "Mean Correlations (Controls)",
  cex.main = 2.25,
  col = "lightsteelblue4",
  xlab = "Vehicle Mean Correlation",
  cex.axis = 1.75,
  cex.lab = 1.75
)
abline(v = mean(meancors$MeanCorr),
       col = "black",
       lwd = 4)
abline(
  v = mean(meancors$MeanCorr) - 3 * sd(meancors$MeanCorr),
  col = "red",
  lwd = 2
)
abline(
  v = mean(meancors$MeanCorr) + 3 * sd(meancors$MeanCorr),
  col = "red",
  lwd = 2
)
dev.off()
### Output Dataset of Correlation Matrix ###
write.csv(meancors, file = "Output Files/RawCountResults/Vehicle_replicate_correlations.csv")

### Examination of Plate Controls Pairwise Correlations ###
pdf(file = "Output Files/RawCountResults/PairwiseCorrs.pdf",
    width = 10,
    height = 8)
par(mar = c(2, 2, 5, 2))
chart.Correlation(log((media_counts[, 13:24] + .5), 2), method = "spearman")
mtext(
  side = 3,
  "Pairwise Correlations of Log2(counts + 0.5) for Vehicle Controls",
  line = 3.75,
  cex = 1.5
)
dev.off()


# Remove Outlier controls and Write cleaned Data Sets for DESeq2 --------------------------------------------------

### WRITE OUT COUNTS HERE FOR DESEQ2 for TOS_S2_DESeq2.R ####
### qc_counts and qc_hash are the cleaned input into DESeq2 normalization and DGE identification that runs next ###
### Outliers are removed based on 3sd from the mean of vehicle control sample correlations ###
toremove <-
  rownames(meancors)[(meancors$MeanCorr < (mean(meancors$MeanCorr) - 3 * (sd(meancors$MeanCorr)))) |
                       (meancors$MeanCorr > (mean(meancors$MeanCorr) + 3 * (sd(meancors$MeanCorr))))]
head(qc_hash)
qc_hash$ID2 <- paste0(qc_hash$Treatment, "_", qc_hash$Replicate)
qc_hash <- subset(qc_hash, ID2 != toremove)
qc_counts <- qc_counts[, rownames(qc_hash)]

if (identical(rownames(qc_hash), colnames(qc_counts))) {
  save(qc_hash, file = "Output Files/Cleaned_data/qc_hash.RData")
  save(qc_counts, file = "Output Files/Cleaned_data/qc_counts.RData")
} else {
  stop("Colnames of count matrix don't match rownames of hash file. Please fix")
}


### Heatmap of Well Counts - Only to be run if you have full plate. Modify otherwise###
### This code used randomly generated numbers to demonstrate #########################

# hash2<-hash
# rownames(hash2)<-paste0("X",rownames(hash2))
# fc <- raw_count_data[,rownames(hash2)]
# hash2$well_reads <- colSums(fc)
# all.letters<-LETTERS[1:26]
# mapds<- hash2 %>% select(Index.Set, X96wp.position,Treatment,Sample.ID,well_reads) %>%
#   mutate(row_letter=substr(X96wp.position,1,1),
#          column_num=substr(X96wp.position,2,3)) 
# mapds$rl <- ifelse(mapds$Index.Set %in% c("A","B"), mapds$row_letter, all.letters[match(mapds$row_letter,all.letters)+8])
# mapds$cn <- ifelse(mapds$Index.Set %in% c("A","C"), mapds$column_num, as.character(as.numeric(mapds$column_num) + 12))
# hm_counts<- dcast(rl~cn,data=mapds,value.var = "well_reads",fill=NA)
# row.names(hm_counts) <- hm_counts[,1]
# trt_matrix<-dcast(rl~cn,data=mapds,value.var="Treatment")    
# sample_matrix<-dcast(rl~cn,data=mapds, value.var="Sample.ID")
# 
# 
# 
# library(pheatmap); library(RColorBrewer)
# b<-as.matrix(hm_counts[,-1])
# b<-ifelse(b<100000,abs(rnorm(384,500000,250000)),b)
# b<-b/100000
# b<-ifelse(b<1,NA,b)
# 
# bcolors<- colorRampPalette(c("yellow","orange","purple"), space = "rgb")(15)   
# 
# pdf(file="heatmap_plate.pdf",height= 4.5, width=10.5)    
# pheatmap((as.matrix(b)), display_numbers = T, fontsize_number=12,number_format="%.1f", cluster_rows=FALSE, cluster_cols = FALSE, color=bcolors, main="Simulated Full Count Matrix (100,000s) \n White Cells < 100,000 Total Reads")
# dev.off()





