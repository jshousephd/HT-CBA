#############################################################################################################
##Author :          John House - BRC;NCSU
##Last Modified :   7.24.17 - John House
##Purpose:          Collect DRC, TCPL with the wide DS and Assess POD Best Fits and Numbers
##Info:             R version MRO 3.4.0
##Input:            complete_ds  LL4_data  tcpl_pod
#############################################################################################################
package.list <-
  c("reshape2", "tidyverse", "RColorBrewer", "stringr")
if (length(setdiff(package.list, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(package.list, rownames(installed.packages())))
}
lapply(package.list, require, character.only = TRUE)
setwd("C:/Users/jshouse/Google_Drive/TempOSeq_Pipeline/ht_cba_manuscript")

load(file = "Output Files/DRM/all_pod.RData")
allpod <- all_pod

### identify the best fit model by AIC
names(allpod)

aic_terms <- grep("AIC", names(allpod), ignore.case = T)
allpod$best.fit <-
  colnames(allpod[aic_terms])[apply(allpod[aic_terms], 1, which.min)]
table(allpod$best.fit)
head(allpod, 10)
### REASSIGN BESTFIT NAME TO MATCH THE NAME OF THE COLUMN CONTAINING THE POD VALUE ###
allpod$best.fit <-
  replace(
    allpod$best.fit,
    which(allpod$best.fit == "aic_ucnst" |
            allpod$best.fit == "aic_dcnst"),
    NA
  )
head(allpod)
allpod$best.fit <-
  replace(allpod$best.fit,
          which(allpod$best.fit == "aic_uhill"),
          "hillup_pod")
head(allpod)
allpod$best.fit <-
  replace(allpod$best.fit,
          which(allpod$best.fit == "aic_ugnls"),
          "gnlsup_pod")
head(allpod)
allpod$best.fit <-
  replace(allpod$best.fit,
          which(allpod$best.fit == "aic_dhill"),
          "hilldown_pod")
head(allpod)
allpod$best.fit <-
  replace(allpod$best.fit,
          which(allpod$best.fit == "aic_dgnls"),
          "gnlsdown_pod")
head(allpod, 10)
allpod$best.fit <-
  replace(allpod$best.fit,
          which(allpod$best.fit == "LL4_AIC"),
          "LL4_POD")
head(allpod, 10)
table(allpod$best.fit)

pod_terms <- grep("POD", names(allpod), ignore.case = T)
allpod$log10pod <-
  as.numeric(allpod[pod_terms][cbind(1:nrow(allpod), match(allpod$best.fit, names(allpod)[pod_terms]))])
head(allpod, 10)
summary(allpod$log10pod)
## place bounds on POD calculation to be the range of the dose response ##
allpod$log10pod <-
  ifelse(allpod$log10pod < -2, -2, allpod$log10pod)
head(allpod)
allpod$log10pod <-
  ifelse(allpod$log10pod > 1, 1, allpod$log10pod)
head(allpod)
### tcpl fits both an UP and Down version of data. Both are same Model type.###
allpod$MinAIC_Model <-
  ifelse(
    allpod$best.fit %in% c("hilldown_pod", "hillup_pod"),
    "3P_Hill",
    ifelse(
      allpod$best.fit %in% c("gnlsdown_pod", "gnlsup_pod"),
      "GainLoss",
      ifelse(allpod$best.fit == "LL4_POD", "4P_Hill", allpod$best.fit)
    )
  )

totalmodels <- sum(!is.na(allpod$MinAIC_Model))
t <- table(allpod$MinAIC_Model)
model_freq <- round(t / totalmodels * 100, 1)
labels <- paste0(model_freq, " %")
pdf("Output Files/DRM/fig6c.pdf",height=5.25,width=8)
par(mar = c(5, 5, 5, 5))
bp <- barplot(
  t,
  ylab = "Winning Model Frequency",
  col = "lightsteelblue4",
  cex.lab = 1.75,
  cex.axis = 1.75,
  cex.main = 2,
  cex.names = 2,
  main = "Dose Response Model Fit Selection",
  ylim = c(0, 1.2 * max(table(
    allpod$MinAIC_Model
  )))
)
text(
  y = t,
  x = bp,
  label = labels,
  pos = 3,
  cex = 2
)
dev.off()
table(allpod$MinAIC_Model)

### Characterize the Dose Respose by Chemical ###
allpod$Treatment <-
  str_sub(allpod$id, start = 1,
          end = (str_locate(allpod$id, "/"))[, 1] - 1)
head(allpod)
allpod$Treatment <- factor(allpod$Treatment)
table(allpod$Treatment)

summarystats <- allpod %>% filter(!is.na(log10pod)) %>%
  group_by(Treatment) %>%
  summarise(
    POD_Values = n(),
    Mean_POD = mean(log10pod, na.rm = T),
    Median_POD = median(log10pod)
  )
pdf("Output Files/DRM/fig6d.pdf",height=5.25,width=8)
par(mar = c(5, 5, 5, 5))
x <- boxplot(
  log10pod ~ Treatment,
  data = allpod,
  main = "Distribution of Calculated POD",
  ylab = "Log10(Dose) uM",
  col = "lightsteelblue4",
  border = "darkorange2",
  cex.names = 2,
  cex = 2,
  cex.lab = 1.75,
  cex.axis = 1.75,
  cex.main = 2
)
means <- summarystats[, 3]
points(means,
       col = "purple",
       pch = 18,
       cex = 2)
dev.off()
