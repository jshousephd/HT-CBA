#############################################################################################################
##Author : John House - BRC;NCSU
##Last Modified : 8/25/17 - John House
##Purpose:    Using tcpl to fit dose response curves and generate POD estimates
##Info:       MRO version 3.4
##Input:      stats dataset from TOS_S4...R
##Output:
##Output:
##NOTE:       https://cran.r-project.org/web/packages/tcpl/vignettes/tcpl_Overview.pdf
#############################################################################################################
package.list <- c("tcpl", "tidyverse", "reshape2", "drc")
if (length(setdiff(package.list, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(package.list, rownames(installed.packages())))
}
lapply(package.list, require, character.only = TRUE)
source(file = "R Scripts/getpod_sim.R")


x1 <- c(-1, 0, 1, rep(-2, 12))
x2 <- c(rep(-1, 2), rep(0, 2), rep(1, 2), rep(-2, 12))
x3 <- c(rep(-1, 3), rep(0, 3), rep(1, 3), rep(-2, 12))
y1 <-
  c(
    rnorm(n = 1, mean = 2, sd = .5),
    rnorm(n = 1, mean = 3, sd = .5),
    rnorm(n = 1, mean = 3.5, sd = .75),
    rnorm(n = 12, mean = 0, sd = .75)
  )
y2 <-
  c(
    rnorm(n = 2, mean = 2, sd = .5),
    rnorm(n = 2, mean = 3, sd = .5),
    rnorm(n = 2, mean = 3.5, sd = .75),
    rnorm(n = 12, mean = 0, sd = .75)
  )
y3 <-
  c(
    rnorm(n = 3, mean = 2, sd = .5),
    rnorm(n = 3, mean = 3, sd = .5),
    rnorm(n = 3, mean = 3.5, sd = .75),
    rnorm(n = 12, mean = 0, sd = .75)
  )

###y3 vector for manuscript demonstration ###
# y3 <- c(1.787, 1.566, 1.916, 3.01, 3.85, 3.586, 3.113, 1.655, 2.862, 0.394, -1.209, 0.356, 0.068, 0.433, 0.906, -0.397, -0.373, -0.243,  0.055, -0.855, -0.72)



mean_y3 <- mean(y3[10:21])
mean_cent_y3 <- y3 - mean_y3
mydf <- data.frame(x = x3, y = mean_cent_y3)

x = x3
y = mean_cent_y3
mysd = sd(y[x == min(x)])
mymu = mean(y[x == min(x)])
target = mymu + mysd
myLL = drm(y ~ x,
           data = mydf,
           logDose = 10,
           fct = LL.4(names = c(
             "Slope", "Lower Limit", "Upper Limit", "ED50"
           )))
myED = try(ED(myLL, target, type = "absolute", interval = "delta"))
LL4_AIC <- AIC(myLL)

LL4_RSS <- modelFit(myLL)[2, 2]
up_params <-
  tcplFit(
    logc = x3,
    resp = mean_cent_y3,
    bmad = .01,
    force.fit = TRUE
  )
vsd <-
  sd(mean_cent_y3[10:21]) ## concentration at which to calculate corresponding dose
hillup_pod <-
  tcplHillConc(
    val = vsd,
    bt = 0,
    tp = up_params$hill_tp,
    ga = up_params$hill_ga,
    gw = up_params$hill_gw
  )
gnlsup_pod <-
  tcplHillConc(
    val = vsd,
    bt = 0,
    tp = up_params$gnls_tp,
    ga = up_params$gnls_ga,
    gw = up_params$gnls_gw
  )

pdf(file="Output Files/DRM/fig6a.pdf",height=5.25,width=8)
par(mar = c(6, 6, 6, 6))
plot(
  x = x3,
  y = mean_cent_y3,
  main = "Concentration Response Fitting (tcpl)",
  col = "red",
  cex = 2,
  xlab = paste0("\n", "Log10(Dose) uM"),
  ylab = "Log2 Control Centered Counts",
  pch = 7,
  cex.axis = 1.75,
  cex.lab = 1.75,
  cex.main = 2
)
box(lwd = 2)
tcplAddModel(
  pars = up_params,
  modl = "hill",
  col = "blue",
  lwd = 2
)
tcplAddModel(
  pars = up_params,
  modl = "cnst",
  col = "black",
  lwd = 2
)
tcplAddModel(
  pars = up_params,
  modl = "gnls",
  col = "red",
  lwd = 2
)

abline(h = vsd,
       col = "purple",
       lwd = 2,
       lty = 2)
abline(v = hillup_pod, col = "blue", lwd = 2)
abline(v = gnlsup_pod, col = "red", lwd = 2)
mtext(paste0(
  "AICC=",
  round(up_params$cnst_aic, 2),
  " ; AICH=",
  round(up_params$hill_aic, 2),
  " ; AICGL=",
  round(up_params$gnls_aic)
),
cex = 1.75)
mtext(
  paste0(
    "3P_Hill POD= ",
    round(hillup_pod, 2),
    " ; ",
    "GainLoss POD = ",
    round(gnlsup_pod, 2)
  ),
  side = 1,
  line = 4.5,
  cex = 1.75,
  col = "red"
)

legend(
  "bottomright",
  horiz = F,
  inset = c(0, 0),
  c("BSD", "Cnst", "3PHill", "Gnls"),
  lwd = 2,
  col = c("purple", "black", "blue", "red"),
  cex = 1.5,
  bty = "n",
  lty = c(1, 1, 1),
  text.col = c("purple", "black", "blue", "red"),
  ncol = 2,
  x.intersp = .3,
  y.intersp = .7
)
dev.off()

pdf(file="Output Files/DRM/fig6b.pdf",height=5.25,width=8)
par(mar = c(6, 6, 6, 6))
plot(
  myLL,
  broken = FALSE,
  type = "all",
  col = 2,
  pch = 7,
  add = F,
  lwd = 2,
  cex = 2,
  xlab = "Log10(Dose) uM",
  ylab = "Log2 Control Centered Counts",
  cex.axis = 1.75,
  cex.lab = 1.75,
  cex.main = 2,
  main = paste0("Concentration Response Fitting 4P_Hill")
)
box(lwd = 2)

abline(h = target,
       col = "purple",
       lwd = 2,
       lty = 2)
abline(v = log10(myED[1]),
       col = "red",
       lwd = 2)
abline(h = 0, col = "black", lwd = 2)
mtext(paste0("4PHill.AIC=", round(LL4_AIC, 2)), cex = 1.75)
mtext(
  paste0("4P_Hill POD = ", round(log10(myED[1]), 2), " uM"),
  side = 1,
  line = 4.5,
  cex = 1.75,
  col = "red"
)
dev.off()




#                 plot(myLL, broken=FALSE, xlab="Dose (mM)", ylab="Log2FC",type="all",col=2,cex=2,pch=7,)
#                 title(main=paste0("LL.4: "),sub=paste0("log10(LL.4_POD) = ",round(log10(myED[1]),2)," mM"))
#                 mtext(paste0("LL4.AIC=",round(LL4_AIC,2),"; ResVar=",round(LL4_RSS,2)),cex=.7)
#                     abline(h=target,col="purple",lwd=2,lty=2)
#                     abline(v=log10(myED[1]),col="red",lwd=2)
#                     abline(h=0,col="black",lwd=2)
#                     abline(h=log10(myED[1]),col="red",lwd=2)
#
# legend("bottom", horiz=T, inset=c(0,0),c("4PHill","Cnst","3PHill","Gnls"),lwd=2, col=c("purple","black","blue","red"), cex=.8,bty="n",lty=c(1,1,1),text.col=c("purple","black","blue","red"))
