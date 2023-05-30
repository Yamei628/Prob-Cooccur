################################################################################
# This R-script provides the functions for conducting the tests described in   #
# Chang et al. (2023).                                                         #
################################################################################


################################################################################
#                                                                              #
#            Chi-squared test of independence                                  #
#                                                                              #
################################################################################
Chi.t=function(a,b){
  ####################################################################################
  # a: vector of 0 and 1 representing absence and presence of Species-A, respectively#
  # b: vector of 0 and 1 representing absence and presence of Species-B, respectively#
  ####################################################################################
  
	n1=sum(a) # Number of cells in which Species-A is present
	n2=sum(b) # Number of cells in which Species-B is present
	N=length(a) # Total number of cells/pixels in the study window
  ##############################################################################
  # Compute Expected number                                                    #
	##############################################################################
  E1=n1*n2/N # Both species A and B are present
  E2=(N-n1)*(N-n2)/N # Neither of them is present
	E3=n1*(N-n2)/N # Species A is present but Species B is absent
  E4=(N-n1)*n2/N # Species B is present but Species A is absent
  EE=c(E1,E2,E3,E4) # Vector of expected numbers
  ##############################################################################
  # Compute Observed number                                                    #
  ##############################################################################
	O1=sum(a*b) # Number of cells in which both species A and B are present
	O2=sum((1-a)*(1-b)) # Number of cells in which neither of them is present
	O3=sum(a*(1-b)) # Number of cells in which A is present but B is absent
	O4=sum((1-a)*b) # Number of cells in which B is present but A is absent
  OO=c(O1,O2,O3,O4) # Vector of Observed counts
  ##############################################################################
  # Compute the Chi-squared test statistic                                     #
  ##############################################################################
    temp.j=sum(EE<10)   
    
    chi.s1 = sum( (OO-EE)^2/EE ) # without Yate's correction

    temp=abs(OO-EE)-0.5
    temp1=temp[temp>0]	
    chi.s2=sum(temp1^2/EE[temp>0]) #with Yate's correction

    chi.s = chi.s1*(temp.j==0)+chi.s2*(temp.j>0)
    p=pchisq(chi.s,1,lower.tail=F)
    return(p)
}


################################################################################
#                                                                              #
#   Binomial and Poisson tests for testing positive (or negative) association  #
#                                                                              #  
################################################################################    
BP.t=function(a, b, posi.test=TRUE){
  ####################################################################################
  # a: vector of 0 and 1 representing absence and presence of Species-A, respectively#
  # b: vector of 0 and 1 representing absence and presence of Species-B, respectively#
  # posi.test: When true, the alternative hypothesis of positive association is      #
  #            tested, and when false, the alternative hypothesis of negative        #
  #            association is tested.                                                #
  ####################################################################################
	n1=sum(a) # Number of cells in which Species-A is present
	n2=sum(b) # Number of cells in which Species-B is present
	N=length(a) # Total number of cells/pixels in the study window
  ##############################################################################
	# Compute Expected number                                                    #
	##############################################################################
  E1=n1*n2/N # Both species A and B are present
  E2=(N-n1)*(N-n2)/N # Neither of them is present
	E3=n1*(N-n2)/N # Species A is present but Species B is absent
  E4=(N-n1)*n2/N # Species B is present but Species A is absent
  ##############################################################################
  # Compute Observed number                                                    #
  ##############################################################################
	O1=sum(a*b) # Number of cells in which both species A and B are present
	O2=sum((1-a)*(1-b)) # Number of cells in which neither of them is present
	O3=sum(a*(1-b)) # Number of cells in which A is present but B is absent
	O4=sum((1-a)*b)  # Number of cells in which B is present but A is absent
  ##############################################################################
	# Compute p-values using the tests based on Binomial distribution            #
	##############################################################################
  p1.p=pbinom(O1-1, N, E1/N, lower.tail = FALSE) #P1 for testing positive association
  p1.n=pbinom(O1, N, E1/N, lower.tail = TRUE)    #P1 for tesing negative association
  p1=p1.p*posi.test+p1.n*(1-posi.test)           #P1

  p22.p=pbinom(O2-1, N, E2/N, lower.tail = FALSE) 
  p22.n=pbinom(O2, N, E2/N, lower.tail = TRUE)    
  p2.p=p1.p*p22.p                           #P2 for testing positive association
  p2.n=p1.n*p22.n                           #P2 for tesing negative association
  p2=p2.p*posi.test+p2.n*(1-posi.test)      #P2 

  p31.p=pbinom(O3 , N, E3/N, lower.tail = TRUE)
  p32.p=pbinom(O4 , N, E4/N, lower.tail = TRUE)
  p3.p=p31.p*p32.p             #P3 for testing positive association
  p31.n=pbinom(O3-1, N, E3/N, lower.tail = FALSE) 
  p32.n=pbinom(O4-1, N, E4/N, lower.tail = FALSE) 
  p3.n=p31.n*p32.n             #P3 for tesing negative association
  p3=p3.p*posi.test+p3.n*(1-posi.test)      #P3 
  ##############################################################################
  # Compute p-values using the tests based on Poisson distribution             #
  ##############################################################################
  p4.p=ppois(O1-1, E1, lower.tail = FALSE) #P4 for testing positive association
  p4.n=ppois(O1, E1, lower.tail = TRUE)    #P4 for tesing negative association
  p4=p4.p*posi.test+p4.n*(1-posi.test)           #P4

  p5.p=ppois(O1+O2-1,E1+E2, lower.tail = FALSE) #P4 for testing positive association
  p5.n=ppois(O1+O2, E1+E2, lower.tail = TRUE)    #P4 for tesing negative association
  p5=p5.p*posi.test+p5.n*(1-posi.test)           #P4

  p6.p=ppois(O3+O4, E3+E4, lower.tail = TRUE) #P4 for testing positive association
  p6.n=ppois(O3+O4-1, E3+E4 , lower.tail =FALSE )    #P4 for tesing negative association
  p6=p6.p*posi.test+p6.n*(1-posi.test)           #P4
  
  p.all=c(p1,p2,p3,p4,p5,p6)
  return(round(p.all,digits=5))
}


################################################################################
#                                                                              #
#   Veech's (2013) probabilistic method                                        #
#                                                                              #  
################################################################################  
Veech.t=function(a,b,posi.test=TRUE){
	N=length(a)
	N1=sum(a)
	N2=sum(b)
	Q_obs=sum(a*b)
	
	j.int=max(0,(N1+N2-N))
	j.end=min(N1,N2)
      
      p_lt=NULL
	p_gt=NULL

	for (j in j.int:(Q_obs-1)){
   	pj=choose(N,j)*choose((N-j),(N2-j))*choose((N-N2),(N1-j))/choose(N,N2)/choose(N,N1)
	p_lt=sum(p_lt,pj)
	}


      for (j in (Q_obs+1):j.end){
   	pj=choose(N,j)*choose((N-j),(N2-j))*choose((N-N2),(N1-j))/choose(N,N2)/choose(N,N1)
	p_gt=sum(p_gt,pj)
	}

	j=Q_obs
        pj=choose(N,j)*choose((N-j),(N2-j))*choose((N-N2),(N1-j))/choose(N,N2)/choose(N,N1)
	
      p_lt=p_lt*(Q_obs>j.int)
      p_et=pj*(Q_obs>=j.int)*(Q_obs<=j.end)
      p_gt=p_gt*(Q_obs<j.end)

      p.n=p_lt+p_et
      p.p=p_gt+p_et 
      p.v=p.p*posi.test+p.n*(1-posi.test) 
	return(p.v)
}

