decompose_plots <- function(rules0, rules, data) {
  x11(width=12, height=6)
  par(mfrow=c(2,4),mar= c(4, 4.5, 3.5, 2) + 0.1)
  k=1
  for(i in rules0){
    decompose_rule(rules[i],data,k)
    k=k+1
  }
}

decompose_rule <- function(rule, data,k){

  r = strsplit(rule,'&')[[1]]
  drules <- c(rule,r)
  n = length(drules)
  if(n==3){
    drdata <- generate_ruledata(data,drules)
  }else if(n==4){

    orig <- drules[1]
    r1 <- paste(drules[3],'&',drules[4])
    r2 <- paste(drules[2],'&',drules[4])
    r3 <- paste(drules[2],'&',drules[3])
    drules1 <- c(orig,r1,r2,r3)
    drdata <- generate_ruledata(data,drules1)
  }

  if(n==3){
    drdata_0 = drdata[drdata[,4]==0,][,c(1,2,4)]
    drdata_1 = drdata[drdata[,4]==1,][,c(1,2,4)]
    drdata_2 = drdata[drdata[,6]==1,][,c(1,2,6)]
    drdata_3 = drdata[drdata[,5]==1,][,c(1,2,5)]
    nanb <- which(drdata[,6] == 0 & drdata[,5] == 0)
    drdata_nanb <- drdata[nanb,]
    makeplot3(drdata_0,drdata_1,drdata_2,drdata_3,k)
    # makeplot3b(drdata_0,drdata_1,drdata_2,drdata_3,drdata_nanb,k)


  }else if(n==4){

    drdata_0 = drdata[drdata[,4]==0,][,c(1,2,4)] # All 3
    drdata_1 = drdata[drdata[,4]==1,][,c(1,2,4)] # Not endorsed
    drdata_2 = drdata[drdata[,5]==1,][,c(1,2,5)]  # B&C
    drdata_3 = drdata[drdata[,6]==1,][,c(1,2,6)]  # A&C
    drdata_4 = drdata[drdata[,7]==1,][,c(1,2,7)]  # A&B
    # drdata_4 = drdata[drdata[,8]==1,][,c(1,2,8)]  # B&C
    # drdata_5 = drdata[drdata[,9]==1,][,c(1,2,9)]  # A&C
    # drdata_6 = drdata[drdata[,10]==1,][,c(1,2,10)]  # A&B

    makeplot4(drdata_0,drdata_1,drdata_2,drdata_3,drdata_4,k)
  }
  else if(n==2){
    drdata <- generate_ruledata(data,drules[1])
    plot(survfit(Surv(time,delta)~drdata[,4],data=drdata), conf.int = FALSE,
         col=c(2,3), lwd=3, lty=1,  ylab = "Survival Probability", xlab='Time',
         ylim=range(0,1), xlim=range(0,5), cex.lab = 1.98, cex.axis = 1.75)
    legend("topright",legend=c("Not Endorsed","Endorsed"),col=c(2,3),
           lty=c(1), lwd=2, cex=1, bty = 'n')
    title(main = paste('Decomposition: Rule', k), cex.main = 1.8)
  }
}

makeplot3<-function(drdata_0,drdata_1,drdata_2,drdata_3,k){
  plot(survfit(Surv(time,delta)~drdata_0[,3],data=drdata_0), conf.int = FALSE,
       col=2, lwd=2, lty=1,ylim=range(0,1), xlim=range(0,5),cex.axis = 1.75,axes=FALSE)
  par(new=TRUE)
  plot(survfit(Surv(time,delta)~drdata_1[,3],data=drdata_1), conf.int = FALSE,
       col=3, lwd=2, lty=1,ylim=range(0,1), xlim=range(0,5),cex.axis = 1.75,axes=FALSE)
  par(new=TRUE)
  plot(survfit(Surv(time,delta)~drdata_2[,3],data=drdata_2), conf.int = FALSE,
       col=c(4), lwd=2, lty=1,ylim=range(0,1), xlim=range(0,5),cex.axis = 1.75,axes=FALSE)
  par(new=TRUE)
  plot(survfit(Surv(time,delta)~drdata_3[,3],data=drdata_3), conf.int = FALSE,
       col="orange", lwd=2, lty=1,  ylab = "Survival Probability", xlab='Time', ylim=range(0,1), xlim=range(0,5),
       cex.lab = 1.98, cex.axis = 1.75)
  legend("topright",legend=c("Not Endorsed","Endorsed", "Remove A","Remove B"),col=c(2,3,4,"orange"),
         lty=1, lwd=2, cex=1, bty = 'n')
  title(main = paste('Decomposition: Rule', k), cex.main = 1.8)
}

makeplot4<-function(drdata_0,drdata_1,drdata_2,drdata_3,drdata_4,k){
  plot(survfit(Surv(time,delta)~drdata_0[,3],data=drdata_0), conf.int = FALSE,
       col=c(2), lwd=2, lty=1,ylim=range(0,1), xlim=range(0,5),cex.axis = 1.75,axes=FALSE)
  par(new=TRUE)
  plot(survfit(Surv(time,delta)~drdata_1[,3],data=drdata_1), conf.int = FALSE,
       col=c(3), lwd=2, lty=1,ylim=range(0,1), xlim=range(0,5),cex.axis = 1.75,axes=FALSE)
  par(new=TRUE)
  plot(survfit(Surv(time,delta)~drdata_2[,3],data=drdata_2), conf.int = FALSE,
       col=4, lwd=2, lty=1,ylim=range(0,1), xlim=range(0,5),cex.axis = 1.75,axes=FALSE)
  par(new=TRUE)
  plot(survfit(Surv(time,delta)~drdata_3[,3],data=drdata_3), conf.int = FALSE,
       col=5, lwd=2, lty=1,ylim=range(0,1), xlim=range(0,5),cex.axis = 1.75,axes=FALSE)
  par(new=TRUE)
  plot(survfit(Surv(time,delta)~drdata_4[,3],data=drdata_4), conf.int = FALSE,
       col=c(6), lwd=2, lty=1,  ylab = "Survival Probability", xlab='Time', ylim=range(0,1),
       xlim=range(0,5), cex.lab = 1.98, cex.axis = 1.75)
  legend("topright",legend=c("Not Endorsed","Endorsed", "Remove A","Remove B","Remove C"),col=c(2:6),
         lty=1,lwd=c(2), cex=1, bty = 'n')
  title(main = paste('Decomposition: Rule', k), cex.main = 1.8)
}

