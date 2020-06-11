# Analyze Hx transfer by double kinetic fit.
# GPL 3.0 
# by Jeremy Deuel <jeremy@deuel.ch>
# 25.01.2020
# Description: This R script analyzes a deconvoluted scanning kinetics file of a Hx transfer experiment
# Input: Deconvoluted scanning kinetics file

# Settings
# set working directory
setwd(getSrcDirectory()[1])
# define the variable dcnvFile to the path where the deconvoluted csv file (generated bz scanning kinetics.R) is to be found.
dcnvFile <- "experiments/example.out.1.hpx1.csv"
dcnv <- read.csv(dcnvFile)

# generate plot

# set up plot
maxTime <- max(dcnv$time)
maxTime <- 100 #override
maxConc <- max(sapply(seq(3,dim(dcnv)[2]),function(x) max(dcnv[,x])))
colorScheme <- c(
	"hpx"="#AC021588",
	"fit"="#00000088",
	"HbMet"="#583A2A88",
	"residue"="#756C6D88",
	"precip"="#ED947F88"
)
plot(0,0,type="l",xlim=c(0,maxTime),ylim=c(0,maxConc),xlab="Time [min]",ylab="Concentration [µM]",bty="n",main=dcnvFile)
sapply(seq(3,dim(dcnv)[2]),function(i) {
	if (colnames(dcnv)[i]=="hpx")
		points(dcnv$time,dcnv[[i]],col=colorScheme[colnames(dcnv)[i]],pch=19,cex=0.25)	
	else
		lines(dcnv$time,dcnv[[i]],col=colorScheme[colnames(dcnv)[i]],lwd=2)	
})

fit_biexp <- function(dta,var,start=c(m1=max(dta[[var]])/2,m2=max(dta[[var]])/2,k2=0.01)) {
	m <-nls(as.formula(paste(var,"~m1*(1-exp(-300*time))+m2*(1-exp(-k2*time))")),dta,start=start)
	co <- coef(m)
	return(m)
}

m <- fit_biexp(dcnv,"hpx")
lines(dcnv$time,predict(m),col="#00000088",lwd=2)

legend(maxTime,maxConc/2,c("Measured Hx","Fitted Hx","Measured MetHb","Fit Residue","Precip"),lwd=c(2,2,2,2,2),lty=c(3,1,1,1,1),col=colorScheme,bty="n",xjust=1,yjust=0.5)



cat("RESULT for",dcnvFile)
cat("======")
cat("Free heme transfer with kbHx 50s-1")
cat("Amount of heme: ",round(coef(m)[1],digits=2),"µM")
cat("Heme release by HP")
cat("Slow constant: kslow=",round(coef(m)[3]/60*1000,digits=2),"10-3/s, half life =",round(log(2)/coef(m)[3]*60,0),"s")
cat("Amount of heme: ",round(coef(m)[2],digits=2),"µM")
cat("Fraction of quick released heme: ",round(100*coef(m)[1]/sum(coef(m)[c(1,2)]),digits=1),"%")
print(
results <- c(
	"kr [E-3 s-1]"=as.numeric(coef(m)[3]/60*1000),
	"Kd [µM]"=as.numeric(coef(m)[1]^2/coef(m)[2]),
	"kb [E-3 s-1 µM-1]"=as.numeric(1000*coef(m)[3]/60/(coef(m)[1]^2/coef(m)[2]))
)
)



#draw these kinetics

t <- 1:maxTime
lines(t,coef(m)[1]*(1-exp(-300*t)),lty=2,lwd=1)
lines(t,coef(m)[1]+coef(m)[2]*(1-exp(-coef(m)[3]*t)),lty=2,lwd=1)