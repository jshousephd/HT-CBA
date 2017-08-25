#############################################################################################################
##Author :        John House - BRC;NCSU
##Last Modified : 6/28/17 - John House
##Purpose:        Using tcpl to fit dose response curves and generate POD estimates
##Info:           R version 3.4
##Input:          stats dataset from TOS_S4...R
##Output:
##Output:
##NOTE:    https://cran.r-project.org/web/packages/tcpl/vignettes/tcpl_Overview.pdf
#############################################################################################################
package.list <- c("tcpl", "tidyverse", "reshape2")
if (length(setdiff(package.list, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(package.list, rownames(installed.packages())))
}
lapply(package.list, require, character.only = TRUE)
setwd("C:/Users/jshouse/Google_Drive/TempOSeq_Pipeline/ht_cba_manuscript")

load(file = "Output Files/Cleaned_data/stats.RData")
head(stats)
stats$id <- paste0(stats$treatment, "/", stats$probe)
## Select data with overall trend ##
adj_pval_cutoff <-
  .15  ####change here if a more stringent Cutoff is desired ####
sp_flag_stats <- filter(stats, mcc_qval_all <= adj_pval_cutoff)
## Select where trt != controls via wilcoxon rank sum test for differences, or if they are equal, and trts are significant select those too ##
wilcox_flag <-
  filter(sp_flag_stats, (Wilcoxon_qval <= adj_pval_cutoff) |
           ((Wilcoxon_qval <= adj_pval_cutoff) &
              (mcc_qval_trt < adj_pval_cutoff)
           ))

## Examine if any cases exist where (overall trend signficant, vehicles != treatments, but trt trend was not significant)
flagged_data_report <-
  filter(sp_flag_stats,
         Wilcoxon_qval > adj_pval_cutoff &
           mcc_qval_trt > adj_pval_cutoff)
save(flagged_data_report, file = "Output Files/DRM/flagged_data_report.RData")
write.csv(flagged_data_report, file = "Output Files/DRM/flagged_data_report.csv")

## CHANGE ALL VALUES OF DOSE AND CONTROLS TO BE = LOG2((COUNT+0.5)/AVG_VEHICLE_COUNTS)###
veh_only_idx <- grep("VEH", names(wilcox_flag), ignore.case = T)
trt_only_idx <- grep("dose", names(wilcox_flag), ignore.case = T)
allcounts_idx <-
  c(trt_only_idx, veh_only_idx, length(names(wilcox_flag)))

### Continue with Dose Response Modeling ###
head(wilcox_flag)
log_comp_ds <- wilcox_flag[, allcounts_idx]
head(log_comp_ds)
log_comp_ds$veh_avg_counts <-
  apply(log_comp_ds[, veh_only_idx], 1, function(x)
    mean(x, na.rm = T))
head(log_comp_ds)

## TCPL requires log2(mean centered counts ###
for (i in c(trt_only_idx, veh_only_idx)) {
  log_comp_ds[i] <-
    apply(log_comp_ds[c(i, length(log_comp_ds))], 1, function(x)
      log((x[1] + .5) / x[2], 2))
}

head(log_comp_ds)

### save(log_comp_ds,file="Output Files/log_comp_ds.RDATA")
aic_ucnst = NULL
aic_uhill = NULL
aic_ugnls = NULL
hillup_pod = NULL
gnlsup_pod = NULL
aic_dcnst = NULL
aic_dhill = NULL
aic_dgnls = NULL
hilldown_pod = NULL
gnlsdown_pod = NULL
res.var_ucnst = NULL
res.var_uhill = NULL
res.var_ugnls = NULL
res.var_dcnst = NULL
res.var_dhill = NULL
res.var_dgnls = NULL
id <- NULL

mynames <- names(log_comp_ds)
dose_vector <-
  sapply(strsplit(mynames, split = "_", fixed = T)[c(trt_only_idx, veh_only_idx)], function(x)
    (as.numeric(x[2])))
ld <- log10(replace(dose_vector, dose_vector == 0, .01))

ptm <- proc.time()
for (i in 1:nrow(log_comp_ds)) {
  if (i %% 100 == 0) {
    print(paste0(i, " of ", nrow(log_comp_ds)))
  }
  mydata <- log_comp_ds[i,]
  response <-
    as.vector(as.matrix(mydata[1, c(trt_only_idx, veh_only_idx)]))
  neg_response = -(response)
  vsd <- sd(mydata[veh_only_idx])
  up_params <-
    try(tcplFit(
      logc = ld,
      resp = response,
      bmad = 0,
      force.fit = TRUE
    ))
  down_params <-
    try(tcplFit(
      logc = ld,
      resp = neg_response,
      bmad = 0,
      force.fit = TRUE
    ))
  hillup_pod[i] <-
    try(tcplHillConc(
      val = vsd,
      bt = 0,
      tp = up_params$hill_tp,
      ga = up_params$hill_ga,
      gw = up_params$hill_gw
    ))
  gnlsup_pod[i] <-
    try(tcplHillConc(
      val = vsd,
      bt = 0,
      tp = up_params$gnls_tp,
      ga = up_params$gnls_ga,
      gw = up_params$gnls_gw
    ))
  hilldown_pod[i] <-
    try(tcplHillConc(
      val = vsd,
      bt = 0,
      tp = down_params$hill_tp,
      ga = down_params$hill_ga,
      gw = down_params$hill_gw
    ))
  gnlsdown_pod[i] <-
    try(tcplHillConc(
      val = vsd,
      bt = 0,
      tp = down_params$gnls_tp,
      ga = down_params$gnls_ga,
      gw = down_params$gnls_gw
    ))
  
  id[i] <- mydata$id
  aic_ucnst[i] = up_params$cnst_aic
  aic_uhill[i] = up_params$hill_aic
  aic_ugnls[i] = up_params$gnls_aic
  res.var_ucnst[i] = (up_params$cnst_rmse) ^ 2
  res.var_uhill[i] = (up_params$hill_rmse) ^ 2
  res.var_ugnls[i] <- (up_params$gnls_rmse) ^ 2
  aic_dcnst[i] = down_params$cnst_aic
  aic_dhill[i] = down_params$hill_aic
  aic_dgnls[i] = down_params$gnls_aic
  res.var_dcnst[i] = (down_params$cnst_rmse) ^ 2
  res.var_dhill[i] = (down_params$hill_rmse) ^ 2
  res.var_dgnls[i] <- (down_params$gnls_rmse) ^ 2
}

tcpl_pod <-
  data.frame(
    id,
    hillup_pod,
    gnlsup_pod,
    hilldown_pod,
    gnlsdown_pod,
    aic_ucnst,
    aic_uhill,
    aic_ugnls,
    aic_dcnst,
    aic_dhill,
    aic_dgnls,
    res.var_ucnst,
    res.var_uhill,
    res.var_ugnls,
    res.var_dcnst,
    res.var_dhill,
    res.var_dgnls
  )
proc.time() - ptm

save(tcpl_pod, file = "Output Files/DRM/tcpl_pod.RData")
write.csv(tcpl_pod, file = "Output Files/DRM/tcpl_pod.csv")

######################################################################################################
### Fitting a 4-parameter hill function to the data using DRC package ###
######################################################################################################
package.list <- c("drc")
if (length(setdiff(package.list, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(package.list, rownames(installed.packages())))
}
lapply(package.list, require, character.only = TRUE)
### Load wrapper function call to drc ###
source(file = "R Scripts/jsh_getpod.R")

head(log_comp_ds)
myds1 <- log_comp_ds

LL4_POD <-
  NULL
LL4_AIC <- NULL
LL4_PODSE <- NULL
LL4_RSS <- NULL
id <- NULL

ptm <- proc.time()
x <- log10(replace(dose_vector, dose_vector == 0, .01))
for (i in 1:nrow(myds1)) {
  if (i %% 25 == 0) {
    print(paste0(i, " of ", nrow(myds1)))
  }
  mydata <- myds1[i,]
  y = as.vector(as.matrix(mydata[1, c(trt_only_idx, veh_only_idx)]))
  temp1 <- try(getpod(x, y, plot.it = F, show.pod = F))
  LL4_POD[i] <-
    temp1[[1]]
  LL4_AIC[i] <-
    temp1[[5]]
  LL4_PODSE[i] <-
    temp1[[2]]
  LL4_RSS[i] <- temp1[[6]]
  id[i] <- mydata$id
}
proc.time() - ptm
LL4_data <- data.frame(id, LL4_POD, LL4_AIC, LL4_RSS, LL4_PODSE)
save(LL4_data, file = "Output Files/DRM/LL4_data.RData")

### Merge with tcpl data ###
all_pod <- left_join(tcpl_pod, LL4_data, by = "id")
save(all_pod, file = "Output Files/DRM/all_pod.RData")
write.csv(all_pod, file = "Output Files/DRM/all_pod.csv")

### single model examination with Graphs #########
# par(mfrow=c(1,2))
# mydata<-filter(myds1,id=="Dofetilide/ACADM_22566")
# x <- log10(replace(dose_vector, dose_vector == 0, .01))
# y = as.vector(as.matrix(mydata[1,c(trt_only_idx,veh_only_idx)]))
# temp1 <- try(getpod(x,y,plot.it=T,show.pod=T))
