library(spatstat)
library(letsR)
source("D://Co-Occ.R")

##########Simulation function for generating two species with negative association##########
SimNeg = function(c1,c2,M){ 
	pv=NULL
	ncell=NULL
      npoint=NULL
	j1 <- function(x,y){ c1-c1*x }
	j2 <- function(x,y){ c2*x }

 	########## Simulation of M times##########
	for (i in 1:M){
	set.seed(i)
	jh1 <- rpoispp(j1)
	jh2 <- rpoispp(j2)

	j1xy <- cbind(jh1$x, jh1$y)
	j2xy <- cbind(jh2$x, jh2$y)

	a=rep('species1',jh1$n)
	b=rep('species2',jh2$n)
	##########PAM##########
	PAM1 <- lets.presab.points(j1xy,a,xmn=0,xmx=1,ymn=0,ymx=1,remove.cells=FALSE,resol=0.05)
	PAM2 <- lets.presab.points(j2xy,b,xmn=0,xmx=1,ymn=0,ymx=1,remove.cells=FALSE,resol=0.05)

	T1=PAM1$Presence_and_Absence_Matrix[,3]
	T2=PAM2$Presence_and_Absence_Matrix[,3]
	########P-values for testing positive association under Binomial and Possison  test###############
      p=BP.t(T1,T2,FALSE)
	p.vee=Veech.t(T1,T2,FALSE)
	pv=rbind(pv,c(p,p.vee))                                #p-value

	n1n2=c(sum(T1),sum(T2))
	ncell=rbind(ncell,n1n2)
    
 	N1N2=c(jh1$n,jh2$n) 
      npoint=rbind(npoint,N1N2)
	}
	output=cbind(npoint,ncell,pv)
	return(output)
}

#######output for table 4###########################################################################################


 c.intensity=c((1:10)*100,(2:5)*1000)
 M=500
 k=length(c.intensity)
 outp=NULL
 for (i in 1:k){
 	c12=c.intensity[i]
 	temp=SimNeg(c12,c12,M)
      temp11=temp[,5:11]<0.05
      temp1=c(round(apply(temp[,1:4],2,mean)),apply(temp11,2,sum)/500)
	outp=rbind(outp,temp1)
}

outp



#######output for table 5###########################################################################################


 c.intensity=c((1:10)*100,(2:5)*1000)
 M=500
 k=length(c.intensity)
 outp=NULL
 for (i in 1:k){
 	c12=c.intensity[i]
 	temp=SimNeg(400,c12,M)
      temp11=temp[,5:11]<0.05
      temp1=c(round(apply(temp[,1:4],2,mean)),apply(temp11,2,sum)/500)
	outp=rbind(outp,temp1)
}

outp