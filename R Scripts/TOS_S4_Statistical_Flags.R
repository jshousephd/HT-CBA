#############################################################################################################
##Author :          John House - BRC;NCSU
##Last Modified :   6/15/2017 - John House
##Purpose:          Take hash and normcounts from DESeq2; creat wide dataset and compute statistics for Flags
##NOTE:             ** mcc used to calculate permutation based p-values for trend ** added 3/3/2017
##Info:             MRO version 3.4
##Input:            Hash file and count matrix
##Output:           /Output Files/Cleaned_data/stats.RData  for input into Dose Response Modeling
#############################################################################################################

package.list <- c("tidyverse", "reshape2", "mcc", "wrapr")
if (length(setdiff(package.list, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(package.list, rownames(installed.packages())))
}
lapply(package.list, require, character.only = TRUE)
setwd("C:/Users/jshouse/Google_Drive/TempOSeq_Pipeline/ht_cba_manuscript")

load("Output Files/Cleaned_data/qc_hash.RData")
qc_hash$header <- rownames(qc_hash)
load("Output Files/Cleaned_data/normcounts.RData")
normcounts <- as.data.frame(normcounts)
normcounts$X <- rownames(normcounts)
## First Create Tall Dataset using tidyr::gather
counts <- gather(normcounts, "X", "NormCounts")
names(counts) <- c("Gene", "header", "NormCounts")
## Merge hash file attributes and NormCounts ##
full_ds <- counts %>% left_join(. , qc_hash, by = "header") %>%
  mutate(
    DOSE = paste0("Dose_", numeric_dose),
    Treatment_Dose_Rep = paste0(Treatment, "_", numeric_dose, "_", Replicate),
    Treatment = as.character(Treatment),
    trt_gene = paste0(Treatment, "_", Gene)
  ) %>%
  select(trt_gene,
         Gene,
         NormCounts,
         Treatment,
         numeric_dose,
         Treatment_Dose_Rep)
head(full_ds)

trtments <- full_ds %>% select(Treatment) %>%
  filter(Treatment != "VEHICLE") %>%
  unique() %>% .[, 1]

final_ds <- data.frame()
for (i in 1:length(trtments)) {
  temp <- full_ds %>%
    filter(Treatment %in% trtments[i] | Treatment == "VEHICLE") %>%
    dcast(Gene ~ Treatment_Dose_Rep, value.var = "NormCounts") %>%
    mutate(treatment = trtments[i],
           trt_gene = paste0(trtments[i], "/", Gene))
  names(temp)[c(2:4)] <- c("dose_0.1_1", "dose_1_1", "dose_10_1")
  final_ds <- rbind(final_ds, temp)
}
head(final_ds)

rownames(final_ds) <- final_ds$trt_gene
stats_ds <- select(final_ds, -Gene, -treatment, -trt_gene)
# Define Dose Vector Based on Names -------------------------------------------------------------------------------
dosevector <-
  sapply(strsplit(names(stats_ds), split = "_", fixed = T), function(x)
    (as.numeric(x[2])))
# Define zero dose as one distance closer to zero from lowest dose ------------------------------------------------
log_dose_vector <- log10(replace(dosevector, dosevector == 0, .01))
log_dose_trt_vector <- log_dose_vector[1:3]

## Calculate Statistics for Dataset ##
mynames <- names(stats_ds)
veh_only_idx <- seq(1:length(mynames))[regexpr("VEH", mynames) > 0]
trt_only_idx <- seq(1:length(mynames))[regexpr("VEH", mynames) < 0]
all_idx <- seq(1:length(mynames))

## Remove Records with < 5 counts average across all columns ##
stats_ds <- stats_ds[rowMeans(stats_ds) > 5, ]

head(stats_ds)

stats_ds$p_rho <-
  apply(stats_ds[, all_idx], 1, function(x)
    cor.test(x, y = log_dose_vector, method = "pearson")$estimate)
stats_ds$s_rho <-
  apply(stats_ds[, all_idx], 1, function(x)
    cor.test(x, y = log_dose_vector, method = "spearman")$estimate)
stats_ds$p_rho_pval <-
  apply(stats_ds[, all_idx], 1, function(x)
    cor.test(x, y = log_dose_vector, method = "pearson")$p.value)
stats_ds$s_rho_pval <-
  apply(stats_ds[, all_idx], 1, function(x)
    cor.test(x, y = log_dose_vector, method = "spearman")$p.value)
stats_ds$Wilcoxon_pval <-
  apply(stats_ds[, all_idx], 1, function(x)
    wilcox.test(x[trt_only_idx], x[veh_only_idx])[[3]])
stats_ds$spear_trt_pval <-
  apply(stats_ds[, trt_only_idx], 1, function(x)
    cor.test(x, y = log_dose_trt_vector, method = "spearman")$p.value)

### mcc pvalues of trend ###
x <- as.matrix(stats_ds[all_idx])
y <- log_dose_vector
stats_ds$mcc_pval_all <-
  getbetap.A(
    getAmoment.list = getAmoment(x, y, z = NULL),
    A = NULL,
    fix.obs = FALSE
  )[[1]]
stats_ds$mcc_pval_trt <-
  getbetap.A(
    getAmoment.list = getAmoment(x[, 1:3], y[1:3], z = NULL),
    A = NULL,
    fix.obs = FALSE
  )[[1]]
save(stats_ds, file = "Output Files/Cleaned_data/stats_ds_full.RData")
## Calculate Adjusted Pvalues for each pvalue by Treatment. i.e. Qvalue for each 2800 genes of each treatment ##

stats_ds1 <- stats_ds
stats_ds1$treatment <-
  sapply(rownames(stats_ds1), function(x)
    unlist(strsplit(x, "/"))[[1]])
stats_ds1$probe <-
  sapply(rownames(stats_ds1), function(x)
    unlist(strsplit(x, "/"))[[2]])
stats_ds1$gene <-
  sapply(stats_ds1$probe, function(x)
    unlist(strsplit(x, "_"))[[1]])
treatment <- unique(stats_ds1$treatment)


qvalue_trt <- function(ds.name, pval.name) {
  let(
    c(x = pval.name),
    qval <- ds.name %>% dplyr::select(probe, treatment, x) %>%
      group_by(treatment) %>%
      mutate(q = p.adjust(x, method = c("BH")))
  )
  return(as.vector(qval$q))
}

stats_ds1$p_rho_qval <- qvalue_trt(stats_ds1, "p_rho_pval")
stats_ds1$s_rho_qval <- qvalue_trt(stats_ds1, "s_rho_pval")
stats_ds1$Wilcoxon_qval <- qvalue_trt(stats_ds1, "Wilcoxon_pval")
stats_ds1$spear_trt_qval <- qvalue_trt(stats_ds1, "spear_trt_pval")
stats_ds1$mcc_qval_all <- qvalue_trt(stats_ds1, "mcc_pval_all")
stats_ds1$mcc_qval_trt <- qvalue_trt(stats_ds1, "mcc_pval_trt")

stats <- stats_ds1

save(stats, file = "Output Files/Cleaned_data/stats.RData")
save(full_ds, file = "Output Files/Cleaned_data/full_ds.RData")
## THE STATS DATASET --> FEED INTO DOSE RESPONSE MODELING ##

# confirming qvalue_trt code ###
# ttt <- stats_ds1 %>% dplyr::filter(treatment == "Nifedipine") %>% dplyr::select(p_rho_qval)
# tt <- stats_ds1 %>% dplyr::filter(treatment == "Nifedipine") %>% dplyr::select(p_rho_pval)
# identical(p.adjust(tt[,1],method=c("BH")),ttt[,1])
