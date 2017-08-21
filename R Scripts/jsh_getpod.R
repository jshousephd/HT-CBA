### WRAPPER FUNCTION TO DRC TO FIT A 4 PARAMETER HILL FUNCTION, CALCULATE POD AND RETURN STAT###

getpod=function(x,y,plot.it=F,show.pod=F){
  pod=NA;  podSE=NA;  podlower=NA;  podupper=NA; LL4_AIC=NA; LL4_RSS=NA
  mydf=data.frame(x,y)
  mysd=sd(y[x==min(x)])
  mymu=mean(y[x==min(x)])
  myLL = try(drm(y ~ x, data = mydf,logDose=10, fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50"))))
  
  #  myMod_sum = try(as.numeric(mselect(myLL,icfct = AIC)))
  if (sum(is.na(myLL) == 0) & regexpr("Error",myLL[1]) == -1){
    mycoef=summary(myLL)$coef[,1]
    #LL4_AIC=myMod_sum[2]
    
    #LL4_resvar=myMod_sum[4]
    LL4_AIC <- AIC(
      myLL);
    LL4_RSS <- modelFit(myLL)[2,2]
    b=mycoef[1]
    c=mycoef[2]
    d=mycoef[3]
    ED50=mycoef[4]
    finf=c*(sign(b)==1)+d*(sign(b)==-1)
    f0=d*(sign(b)==1)+c*(sign(b)==-1)
    if (finf>f0){target=mymu+mysd}
    if (finf<f0){target=mymu-mysd}
    myED=try(ED(myLL,target,type="absolute",interval="delta"))
    pod <- ifelse(myED[1]<.01,NA,ifelse(myED[1]>1,NA,myED[1]))
    #pod=myED[1]
    podSE=myED[2]
    podlower=myED[3]
    podupper=myED[4]
    myED[1]<-round(pod,5)
    
    if (plot.it==T){
      plot(myLL, broken=FALSE, xlab="Dose (mM)", ylab="Log2 (Mean Centered Normalized Counts)",type="all",col=2,cex=1.5,pch=7)
      title(main=paste0("LL.4: ", mydata$id[1]),sub=paste0("log10(LL.4_POD) = ",round(log10(myED[1]),2)," mM"))
      mtext(paste0("LL4.AIC=",round(LL4_AIC,2),"; ResVar=",round(LL4_RSS,2)),cex=.7)
      if (show.pod==T){
        abline(h=target,col="purple",lwd=2,lty=2)
        abline(v=log10(pod),col="red",lwd=2)
        abline(h=0,col="black",lwd=2)
      }
    }
  }
  return(list(pod=log10(pod),podSE=podSE,podlower=podlower,podupper=podupper,LL4_AIC=LL4_AIC,LL4_RSS=LL4_RSS))
  i=i+1
}