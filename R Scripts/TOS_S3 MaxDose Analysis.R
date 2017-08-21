#############################################################################################################
##Author :          John S. House - BRC;NCSU
##Last Modified :   8.16.17 - John House 
##Purpose:          Take MaxDose L2FC Conduct Cluster Analyses
##Info:             MRO version 3.4
##Input:            combined_maxdose from TOS_S2....R
##Output:           ABS_LOG2FC_Boxplot.pdf  Corrmatrix_l2fc.pdf  heatmap_l2fc.pdf  PCA_L2fc.pdf
##Output:
##NOTE:             The visualizations here are on Log2Fold change from DESeq2 Analysis
#############################################################################################################
package.list <-
  c(
    "tidyverse",
    "matrixStats",
    "dendextend",
    "reshape2",
    "corrplot",
    "gplots",
    "d3heatmap",
    "scatterplot3d",
    "plot3D",
    "scatterD3"
  )
if (length(setdiff(package.list, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(package.list, rownames(installed.packages())))
}
lapply(package.list, require, character.only = TRUE)
setwd("C:/Users/jshouse/Google_Drive/TempOSeq_Pipeline/ht_cba_manuscript")
load(file = "Output Files/MaxDoseContrasts/combined_maxdose.RData")
combined_maxdose$absl2fc <- abs(combined_maxdose$log2FoldChange)
### Create Log2FC matrix from combined_maxdose ###
l2fc_m <-
  dcast(combined_maxdose, gene ~ Treatment, value.var = "log2FoldChange")
rownames(l2fc_m) <- l2fc_m[, 1]
l2fc_m <- l2fc_m[, -1]

### Dendrograms ###
pdf(file = "Output Files/MaxDoseContrasts/ABS_LOG2FC_Clustering_Dend.pdf",
    width = 10.5,
    height = 8.0)
col.clus <- hclust(dist(t(l2fc_m)), method = "average")
dend <-
  col.clus %>% as.dendrogram() %>% color_branches(k = 3) %>% hang.dendrogram(hang_height = 5) %>%
  set("branches_lwd", 4) %>% plot(
    main = paste("Average Clustering Log2FC Cardiomyocytes"),
    horiz = TRUE,
    nodePar = list(cex = .007)
  )
dev.off()

### Boxplots ###
pdf(file = "Output Files/MaxDoseContrasts/ABS_LOG2FC_Boxplot.pdf",
    width = 6.3,
    height = 4.84)
par(mar = c(8, 5, 5, 2))
boxplot(
  abs(log2FoldChange) ~ Treatment,
  data = combined_maxdose,
  main = "Expression Changes by Treatment",
  ylab = "Absolute (Log2FC)",
  cex.lab = 1.5,
  cex.main = 1.5,
  cex.axis = 1.25,
  cex = 1.5,
  col = "lightsteelblue4",
  border = "darkorange2",
  las=2
)
dev.off()

#### CORRELATION MATRIX #### ##  https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html
colorscale <-
  colorRampPalette(c("orange", "white", "purple"), space = "rgb")(100)

pdf(file = "Output Files/MaxDoseContrasts/Corrmatrix_l2fc.pdf",
    width = 10.5,
    height = 8.0)
corrplot(
  cor(l2fc_m),
  order = "hclust",
  hclust.method = "average",
  title = "Correlation Matrix of Log2FC at Max Dose iPSC Cardiomyocytes",
  mar = c(2, 2, 2, 2),
  tl.col = "black",
  col = colorscale,
  cl.pos = "b",
  tl.pos = "d",
  method = "number"
)
dev.off()

### Heatmap      ##########################
row.clus <- hclust(dist(l2fc_m), method = "average")
col.clus <- hclust(dist(t(l2fc_m)), method = "average")
colorscale <-
  colorRampPalette(c("purple", "white", "orange"), space = "rgb")(20)
pdf(file = "Output Files/MaxDoseContrasts/heatmap_l2fc.pdf",
    width = 6.3,
    height = 4.84)
heatmap.2(
  as.matrix(l2fc_m),
  Colv = as.dendrogram(col.clus),
  Rowv = as.dendrogram(row.clus),
  col = colorscale,
  density.info = "none",
  labRow = FALSE,
  trace = "none",
  dendrogram = "col",
  srtCol = 90,
  keysize = 1,
  key.title = NA,
  cexCol = 1.5,
  key.xlab = NULL,
  main = "Treatment Clustering of L2FC",
  cex=1,
  margins = c(9, 1)
)
dev.off()

library(scatterplot3d) 
### 3D scatterplot of PCA components for Manuscript
p = prcomp(l2fc_m, center = T)
pc <- as.data.frame(unclass(p$rotation)[, 1:3])
p_var = (round(100 * (p$sdev ^ 2 / sum(p$sdev ^ 2)), 1))[1:3]

pdf(file = "Output Files/MaxDoseContrasts/PCA_L2fc1.pdf",
    width = 6.3,
    height = 4.84)
s3d <- scatterplot3d(
  x = pc$PC1,
  y = pc$PC2,
  z = pc$PC3,
  color = "purple",
  pch=16,
  cex.symbols=2,
  cex.lab = 1.25,
  xlab = paste0("PC1 - ", p_var[1], "% Variance"),
  ylab = paste0("PC2 - ", p_var[2], "% Variance"),
  zlab = paste0("PC3 - ", p_var[3], "% Variance"),
  main = "Treatment PCA (L2FC)", cex.main=1.5
)
text(s3d$xyz.convert(pc$PC1, pc$PC2, pc$PC3), labels = rownames(pc), pos = 4 )
dev.off()

### Summary of how Many Max Dose Responsive Genes Found by Treatment ####
load(file = "Output Files/MaxDoseContrasts/combined_maxdose.RData")
q_threshold = 0.10

q_sig <- combined_maxdose %>%
  dplyr::filter(padj < q_threshold) %>%
  mutate(Direction = ifelse(log2FoldChange < 0, "Up", "Down")) %>%
  group_by(Treatment, Direction) %>%
  summarise(
    Gene_number = n(),
    mean_abslog2fc = mean(abs(log2FoldChange)),
    median_abslog2fc = median(abs(log2FoldChange))
  ) %>%
  mutate(Significance_Level = q_threshold) %>%
  arrange(-Gene_number)
q_sig

total_s <- q_sig %>%
  group_by(Treatment) %>%
  summarise(tot=sum(Gene_number))

sbc <-
  ggplot(data = q_sig, aes(x = Treatment, y = Gene_number)) +
  geom_bar(stat = "identity", aes(fill=Direction)) +
  geom_text(data= q_sig, aes(y=Gene_number, label = Gene_number ),
            size=4, position = position_stack(vjust=.5)) +
  geom_text(data=total_s,
            aes(y=tot, label = tot), size = 4, vjust=-.5) +
  theme_bw()  +
  theme(text = element_text(size=20),
        axis.text.x = element_text(size = 13, angle = 90, hjust = 1, face = "bold")) +
  scale_fill_manual(
    values = c("snow3", "lightsteelblue1"),
    name = "Direction",
    labels = c("Up", "Down")
  ) +
  ggtitle(paste0("Significant DEGs by Treatment")) + 
  theme(plot.title=element_text(face="bold", size = 20)) +
  labs(y = paste0("Differentially Expressed Genes"), x = "") +
  theme(axis.title=element_text(face="bold", size = 16)) +
  expand_limits(y=c(0, 1.05*max(total_s$tot) ))
print(sbc)

pdf(file = "Output Files/MaxDoseContrasts/Hist_MaxDoseDEGs.pdf",
    width = 6.3,
    height = 4.84)
print(sbc)
dev.off()

### Report of Number of DEGs by Treatment at Max Dose ###
write.csv(q_sig, row.names = F, file = "Output Files/MaxDoseContrasts/Significant_Genes_by_Treatment.csv")



