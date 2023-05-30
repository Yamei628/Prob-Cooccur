setwd("d://")
library(letsR)
library(spatstat)
library(cooccur)
source("Co-Occ.R")

bci.full=read.table("bci2015.txt",header=TRUE,sep="\t")
bci.alive=bci.full[bci.full$Status=="alive",]
bci.u=bci.alive[,c("Latin","PX","PY")]
bci.u=na.omit(bci.u)

w<-500   #width of the study area
l<-1000  #length of the study area
g<-5     #width of each cell
r<-g^2   #area of each cell
kkx=table(bci.u$Latin)
ten.d=kkx[rev(order(kkx))[1:8]]
ten.name=names(ten.d) #names of the first 8 abundant species

XX=NULL
for (i in 1:8){
test.d=bci.u[bci.u[,"Latin"]==ten.name[i],]
  species <-test.d[,"Latin"]
  xy <- test.d[,c("PX","PY")]
  PAM<- lets.presab.points(xy, species, xmn =0, xmx =l,ymn =0, ymx =w,remove.cells=FALSE,resol=5)
  T=PAM$Presence_and_Absence_Matrix[,3]
  XX=rbind(XX,T) #PAM of the first ten abundant species
}
dim(XX)

N <- 200*100 #number of the cells


#########Table ########
data.frame(ten.d,present.cells=apply(XX,1,sum),present.rate=apply(XX,1,sum)/N)


#########Table #########
oupt=NULL
for (i in 1:7){
 for (j in (i+1):8){
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
oupt

#########Table #########
oupt=NULL
for (i in 1:7){
 for (j in (i+1):8){
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


for (i in c(1,6,8)){
	test.d=bci.u[bci.u[,"Latin"]==ten.name[i],]
	 species <-test.d[,"Latin"]
	 xy <- test.d[,c("PX","PY")]
      PAM<- lets.presab.points(xy, species, xmn =0, xmx =l,ymn =0, ymx =w,remove.cells=FALSE,resol=g)
	#x11()
	par(mfrow=c(1,1),mai=c(0.5,0.5,0.5,0.5))
	plot(xy,xlab="",ylab="",main="")
	#x11()
	#par(mfrow=c(1,1),mai=c(0.5,0.5,0.5,0.5))	
	plot(PAM,world=FALSE,legend=FALSE)
}