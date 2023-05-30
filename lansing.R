setwd("d://")
library(letsR)
library(spatstat)
library(cooccur)

source("Co-Occ.R")

data(lansing)
data.t=split(lansing)
data.names=names(data.t)
i.end=length(data.names)
w<-1   #width of the study area
l<-1  #length of the study area
g<-0.05     #width of each cell
r<-g^2   #area of each cell

XX=NULL
N.all=NULL
for (i in 1:i.end){
	test.d=lansing[lansing$marks==data.names[i],]
	species=test.d$marks
	xy=cbind(test.d$x,test.d$y)
      PAM<- lets.presab.points(xy, species, xmn =0, xmx =l,ymn =0, ymx =w,remove.cells=FALSE,resol=g)
  	T=PAM$Presence_and_Absence_Matrix[,3]
  	XX=rbind(XX,T) #PAM of the first ten abundant species
	N.all=rbind(N.all,test.d$n)
}

N <- dim(XX)[2] #number of the cells




#########Table #########
data.frame(names=data.names,N.all,present.cells=apply(XX,1,sum),present.rate=apply(XX,1,sum)/N)


#########Table #########
oupt=NULL
i.end=dim(XX)[1]-1
j.end=dim(XX)[1]
for (i in 1:i.end){
 for (j in (i+1):j.end){
      T1=XX[i,]
	T2=XX[j,]
	p.chi=Chi.t(T1,T2)
        p.bp=BP.t(T1,T2)
	temp=c(i,j,p.chi,p.bp)
	oupt=rbind(oupt,temp)
	}
}
cooccur.lansing <- cooccur(mat = XX, type = "spp_site", thresh = TRUE, spp_names = TRUE)
p.veech=prob.table(cooccur.lansing)[,9]
oupt=cbind(oupt,p.veech)


#########Table #########
oupt=NULL
i.end=dim(XX)[1]-1
j.end=dim(XX)[1]
for (i in 1:i.end){
 for (j in (i+1):j.end){
      	T1=XX[i,]
	T2=XX[j,]
	p.chi=Chi.t(T1,T2)
      	p.bp=BP.t(T1,T2,FALSE)
      	temp=c(i,j,p.chi,p.bp)
	oupt=rbind(oupt,temp)
	}
}

cooccur.lansing <- cooccur(mat = XX, type = "spp_site", thresh = TRUE, spp_names = TRUE)
p.veech=prob.table(cooccur.lansing)[,8]
oupt=cbind(oupt,p.veech)
oupt



##########plot ##################


for (i in 1:i.end){
	test.d=lansing[lansing$marks==data.names[i],]
	species=test.d$marks
	xy=cbind(test.d$x,test.d$y)
      PAM<- lets.presab.points(xy, species, xmn =0, xmx =l,ymn =0, ymx =w,remove.cells=FALSE,resol=g)
	x11()
	par(mfrow=c(1,1),mai=c(0,0,0,0))
	plot(unmark(test.d),main="")
	x11()
	par(mfrow=c(1,1),mai=c(0.5,0.5,0.5,0))	
	plot(PAM,y=c(0,1),xlab="",ylab="",legend=FALSE)
}