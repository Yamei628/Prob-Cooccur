library(spatstat)
library(letsR)
source("D://Co-Occ.R")

##########Random dist.##########
SimRan=function(c,M){
	pv=NULL
	ncell=NULL
	########## Simulation of M times##########
	for (i in 1:M){
	set.seed(i)
	jh1 <- runifpoint(c)
	jh2 <- runifpoint(c)

	j1xy <- cbind(jh1$x, jh1$y)
	j2xy <- cbind(jh2$x, jh2$y)

	a=rep('species1',jh1$n)
	b=rep('species2',jh2$n)

	##########PAM##########
	PAM1 <- lets.presab.points(j1xy,a,xmn=0,xmx=1,ymn=0,ymx=1,remove.cells=FALSE,resol=0.05)
	PAM2 <- lets.presab.points(j2xy,b,xmn=0,xmx=1,ymn=0,ymx=1,remove.cells=FALSE,resol=0.05)

	T1=PAM1$Presence_and_Absence_Matrix[,3]
	T2=PAM2$Presence_and_Absence_Matrix[,3]



	########Chi-sq test###############
	p=Chi.t(T1,T2)
	pv=c(pv,p)                               #p-value

	n1n2=c(sum(T1),sum(T2))
	ncell=rbind(ncell,n1n2)

	}
	output=cbind(ncell,pv)
	return(output)
}


#########################output for table 1#################################

 n.points=c((1:5)*10,(1:10)*100)
 M=500
 k=length(n.points)
 outp=NULL
 for (i in 1:k){
 	n.p=n.points[i]
 	temp=SimRan(n.p,M)
      temp1=c(round(apply(temp[,1:2],2,mean)),sum(temp[,3]<=0.05)/500)
	outp=rbind(outp,temp1)
}

outp