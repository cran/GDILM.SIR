#' This function can use to estimate the model parameters using the initial values.
#'
#' @param ITER Number of iterations
#' @param zz Number of Regions
#' @param lambda0 initial value for Spatial dependence
#' @param sigma0 initial value for the precision of spatial random effects
#' @param Di Euclidean distance between susceptible individual and infectious individual
#' @param g Number of rows in the lattice
#' @param nSample Number of individuals in each cell
#' @param d infectious time units
#' @param n total number of individuals
#' @param time time
#' @param tau the set of infectious individuals at time t in the zth area
#' @param lambda a vector containing the length of infectious period
#' @param alpha0 initial value for the intercept
#' @param beta10 initial value for the parameter corresponding to the covariate associated with susceptible individual
#' @param beta20 initial value for the parameter corresponding to the area-level covariates corresponding to area
#' @param cov1 a vector of covariates associated with susceptible individual
#' @param cov2 a vector of area-level covariates corresponding to area
#' @param phi Spatial random effects
#' @param delta0 Spatial parameter
#' @param Nlabel Label for each sample from the area
#' @param D matrix reflecting neighborhood structure
#' @param I Identity matrix
#'
#' @return the estimated values for the model parameters
#' @export
#'
#'
#'
#' @importFrom MASS mvrnorm
#' @importFrom Matrix bdiag
#' @importFrom stats qnorm
#' @importFrom stats runif
#' @importFrom numDeriv hessian
#'
#' @examples Simulation_Finalmodel(2,4,0.2,0.5,
#' matrix(runif(1600,min=4,max=20),nrow=40,byrow=TRUE),2,10,3,40,10,
#' sample(c(0,1),replace=TRUE,size=40),rep(3,40),0.4,1,1,runif(40,0,1),
#' runif(4,0,1),runif(4,min=0,max=1),2,rep(1:4,each=10),
#' matrix(c(0,-1,0,-1,-1,0,-1,-1,0,-1,0,-1,-1,-1,-1,0),nrow=4,byrow=TRUE),
#' diag(4))
#'
#'
Simulation_Finalmodel <- function(ITER,zz,lambda0,sigma0,Di,g,nSample,d,n,time,
                                  tau,lambda,alpha0,beta10,beta20,cov1,cov2,phi,
                                  delta0,Nlabel,D,I){

GePHI=matrix(0,ITER,zz)
SimPAR=matrix(0,nSample,6)
ComInf=array(0,c(6,6,nSample))
SimRate=array(0,c(n,time,zz,nSample))
SimRate1=array(0,c(n,time,zz,nSample))
Coverage=matrix(0,nSample,6)
PHI=matrix(0,nSample,12)
s1=1
SS=c()
MAX=c()

Sim_Information=function(Nlabel,phi,Di,alpha1,beta1,beta2,delta,sigma1,lambda1,time,n,zz,tau,D,cov1,cov2,lambda,I){
  ######################################## alpha #######################################
  Sim_Exp=function(Nlabel,phi,Di,alpha1,delta,lambda,i,zzz,t,beta1,beta2,n,tau,cov1,cov2){
    SumH1=rep(0,n)
    for(j in 1:n){
      if (Nlabel[j]==zzz){
        if(tau[j]<=t & (tau[j]+lambda[j])>t & tau[j]!=0){
          SumH1[j]=Di[i,j]^(-delta)
        }}}
    SumH1=replace(SumH1,is.infinite(SumH1),0)
    SumH2=sum(SumH1)
    SumH3=exp(alpha1+beta1*cov1[i]+beta2*cov2[zzz])*SumH2
    return(SumH3)
  }


  Sim_ExpDeriv1=function(Nlabel,phi,Di,alpha1,delta,lambda,i,zzz,t,beta1,beta2,n,tau,cov1,cov2){
    SumH4=rep(0,n)
    for(j in 1:n){
      if (Nlabel[j]==zzz){
        if(tau[j]<=t & (tau[j]+lambda[j])>t & tau[j]!=0){
          SumH4[j]=Di[i,j]^(-delta)*(log(Di[i,j]))
        }}}
    SumH4=replace(SumH4,is.infinite(SumH4),0)
    SumH5=sum(SumH4)
    SumH6=exp(alpha1+beta1*cov1[i]+beta2*cov2[zzz])*SumH5
    return(SumH6)
  }

  Sim_ExpDeriv2=function(Nlabel,phi,Di,alpha1,delta,lambda,i,zzz,t,beta1,beta2,n,tau,cov1,cov2){
    SumH7=rep(0,n)
    for(j in 1:n){
      if (Nlabel[j]==zzz){
        if(tau[j]<=t & (tau[j]+lambda[j])>t & tau[j]!=0){
          SumH7[j]=Di[i,j]^(-delta)*(log(Di[i,j]))^2
        }}}
    SumH7=replace(SumH7,is.infinite(SumH7),0)
    SumH8=sum(SumH7)
    SumH9=exp(alpha1+beta1*cov1[i]+beta2*cov2[zzz])*SumH8
    return(SumH9)
  }

  Sim_prob=function(Nlabel,phi,Di,alpha1,delta,i,zzz,t,beta1,beta2,n,tau,lambda,cov1,cov2){
    dx1=rep(0,n)
    for(j in 1:n){
      if (Nlabel[j]==zzz){
        if(tau[j]<=t & (tau[j]+lambda[j])>t & tau[j]!=0){
          dx1[j]=Di[i,j]^(-delta)
        }}}
    dx1=replace(dx1,is.infinite(dx1),0)
    dx=sum(dx1)
    prob1=1-exp(-exp(alpha1+beta1*cov1[i]+beta2*cov2[zzz]+phi[zzz])*dx) #P(i,z,t)
    return(prob1)
  }

  Sim_deriv1=function(Nlabel,phi,Di,alpha1,delta,i,zzz,t,beta1,beta2,n,tau,lambda,cov1,cov2){
    (1- Sim_prob(Nlabel,phi,Di,alpha1,delta,i,zzz,t,beta1,beta2,n,tau,lambda,cov1,cov2))/
      Sim_prob(Nlabel,phi,Di,alpha1,delta,i,zzz,t,beta1,beta2,n,tau,lambda,cov1,cov2)*exp(phi[zzz])
  }

  Sim_deriv2=function(Nlabel,phi,Di,alpha1,delta,i,zzz,t,beta1,beta2,n,tau,lambda,cov1,cov2){
    (1- Sim_prob(Nlabel,phi,Di,alpha1,delta,i,zzz,t,beta1,beta2,n,tau,lambda,cov1,cov2))/
      Sim_prob(Nlabel,phi,Di,alpha1,delta,i,zzz,t,beta1,beta2,n,tau,lambda,cov1,cov2)^2*exp(2*phi[zzz])
  }

  Sim_Estlambda=function(phi,D,I){
    as.numeric(t(phi)%*%(I-D)%*%phi)
  }

  Sim_Estsigma=function(phi,D,lambda1,I){
    as.numeric(t(phi)%*%(lambda1*D+(1-lambda1)*I)%*%phi)
  }

  #######################################################################################################3

  A1=rep(0,time)
  for(t in 1:time){
    A2=rep(0,n)
    for(i in 1:n){
      for(zzz in 1:zz){
        if (Nlabel[i]==zzz){
          if (tau[i]>(t+1)|tau[i]== 0){
            A2[i]=-as.numeric(Sim_Exp(Nlabel,phi,Di,alpha1,delta,lambda,i,zzz,t,beta1,beta2,n,tau,cov1,cov2)*exp(phi[zzz]))
          }}}}
    A1[t]=sum(A2)
  }
  SusA1=sum(A1)
  ###################### infectious period ###################
  A3=rep(0,time)
  for(t in 1:time){
    A4=rep(0,n)
    for(i in 1:n){
      for(zzz in 1:zz){
        if(Nlabel[i]==zzz){
          if(tau[i]==(t+1)){
            A4[i]=as.numeric(Sim_Exp(Nlabel,phi,Di,alpha1,delta,lambda,i,zzz,t,beta1,beta2,n,tau,cov1,cov2)*Sim_deriv1(Nlabel,phi,Di,alpha1,delta,i,zzz,t,beta1,beta2,n,tau,lambda,cov1,cov2))
          }
        }}}
    A3[t]=sum(A4)
  }
  InfA3=sum(A3,na.rm=T)
  EndA3=SusA1+InfA3
  ######################### Second Der...########################
  ##################### infectious period ###################
  A5=rep(0,time)
  for(t in 1:time){
    A6=rep(0,n)
    for(i in 1:n){
      for(zzz in 1:zz){
        if (Nlabel[i]==zzz){
          if(tau[i]==(t+1)){
            A6[i]=-as.numeric((Sim_Exp(Nlabel,phi,Di,alpha1,delta,lambda,i,zzz,t,beta1,beta2,n,tau,cov1,cov2))^2*Sim_deriv2(Nlabel,phi,Di,alpha1,delta,i,zzz,t,beta1,beta2,n,tau,lambda,cov1,cov2))
          }}}}
    A5[t]=sum(A6)
  }
  InfA5=sum(A5,na.rm=T)
  EndA5=EndA3+InfA5
  InfAlpha=EndA5
  #####################################beta ########################
  ##################### infectious period ###################
  BB1=rep(0,time)
  for(t in 1:time){
    BB2=rep(0,n)
    for(i in 1:n){
      for(zzz in 1:zz){
        if (Nlabel[i]==zzz){
          if (tau[i]>(t+1)|tau[i]== 0){
            BB2[i]=-
              as.numeric(cov1[i]^2*Sim_Exp(Nlabel,phi,Di,alpha1,delta,lambda,i,zzz,t,beta1,beta2,n,tau,cov1,cov2)*exp(phi[zzz]))
          }}}}
    BB1[t]=sum(BB2)
  }
  SusBB1=sum(BB1)
  ###################### infectious period ###################
  BB3=rep(0,time)
  for(t in 1:time){
    BB4=rep(0,n)
    for(i in 1:n){
      for(zzz in 1:zz){
        if(Nlabel[i]==zzz){
          if(tau[i]==(t+1)){
            BB4[i]=as.numeric(cov1[i]^2*Sim_Exp(Nlabel,phi,Di,alpha1,delta,lambda,
                                                i,zzz,t,beta1,beta2,n,tau,cov1,cov2)*Sim_deriv1(Nlabel,phi,Di,alpha1,delta,i,zzz,
                                                                                                t,beta1,beta2,n,tau,lambda,cov1,cov2))
          }}}}
    BB3[t]=sum(BB4)
  }
  InfBB3=sum(BB3,na.rm=T)
  EndBB3=SusBB1+InfBB3
  B5=rep(0,time)
  for(t in 1:time){
    B6=rep(0,n)
    for(i in 1:n){
      for(zzz in 1:zz){
        if (Nlabel[i]==zzz){
          if(tau[i]==(t+1)){
            B6[i]=-
              as.numeric(cov1[i]^2*(Sim_Exp(Nlabel,phi,Di,alpha1,delta,lambda,i,zzz,t,beta1,beta2,n,tau,cov1,cov2))^2*Sim_deriv2(Nlabel,phi,Di,alpha1,delta,i,zzz,t,beta1,beta2,n,tau,lambda,cov1,cov2))
          }}}}
    B5[t]=sum(B6)
  }
  InfB5=sum(B5,na.rm=T)
  InfBeta1=EndBB3+InfB5
  #####################################beta2 ########################
  ##################### infectious period ###################
  YY1=rep(0,time)
  for(t in 1:time){
    YY2=rep(0,n)
    for(i in 1:n){
      for(zzz in 1:zz){
        if (Nlabel[i]==zzz){
          if (tau[i]>(t+1)|tau[i]== 0){
            YY2[i]=-
              as.numeric(cov2[zzz]^2*Sim_Exp(Nlabel,phi,Di,alpha1,delta,lambda,i,zzz
                                             ,t,beta1,beta2,n,tau,cov1,cov2)*exp(phi[zzz]))
          }}}}
    YY1[t]=sum(YY2)
  }
  SusYY1=sum(YY1)
  ###################### infectious period ###################
  YY3=rep(0,time)
  for(t in 1:time){
    YY4=rep(0,n)
    for(i in 1:n){
      for(zzz in 1:zz){
        if(Nlabel[i]==zzz){
          if(tau[i]==(t+1)){
            YY4[i]=as.numeric(cov2[zzz]^2*Sim_Exp(Nlabel,phi,Di,alpha1,delta,lambda,i,zzz,t,beta1,beta2,n,tau,cov1,cov2)*Sim_deriv1(Nlabel,phi,Di,alpha1,delta,i,zzz,t,beta1,beta2,n,tau,lambda,cov1,cov2))
          }}}}
    YY3[t]=sum(YY4)
  }
  InfYY3=sum(YY3,na.rm=T)
  EndYY3=SusYY1+InfYY3
  Y5=rep(0,time)
  for(t in 1:time){
    Y6=rep(0,n)
    for(i in 1:n){
      for(zzz in 1:zz){
        if (Nlabel[i]==zzz){
          if(tau[i]==(t+1)){
            Y6[i]=-
              as.numeric(cov2[zzz]^2*(Sim_Exp(Nlabel,phi,Di,alpha1,delta,lambda,i,zzz,t,beta1,beta2,n,tau,cov1,cov2))^2*Sim_deriv2(Nlabel,phi,Di,alpha1,delta,i,zzz,t
                                                                                                                                   ,beta1,beta2,n,tau,lambda,cov1,cov2))
          }}}}
    Y5[t]=sum(Y6)
  }
  InfY5=sum(Y5,na.rm=T)
  InfBeta2=EndYY3+InfY5
  #####################################beta1 beta2########################
  ##################### infectious period ###################
  KK1=rep(0,time)
  for(t in 1:time){
    KK2=rep(0,n)
    for(i in 1:n){
      for(zzz in 1:zz){
        if (Nlabel[i]==zzz){
          if (tau[i]>(t+1)|tau[i]== 0){
            KK2[i]=-
              as.numeric(cov1[i]*cov2[zzz]*Sim_Exp(Nlabel,phi,Di,alpha1,delta,lambda
                                                   ,i,zzz,t,beta1,beta2,n,tau,cov1,cov2)*exp(phi[zzz]))
          }}}}
    KK1[t]=sum(KK2)
  }
  SusKK1=sum(KK1)
  ###################### infectious period ###################
  KK3=rep(0,time)
  for(t in 1:time){
    KK4=rep(0,n)
    for(i in 1:n){
      for(zzz in 1:zz){
        if(Nlabel[i]==zzz){
          if(tau[i]==(t+1)){
            KK4[i]=as.numeric(cov1[i]*cov2[zzz]*Sim_Exp(Nlabel,phi,Di,alpha1,delta,lambda,i,zzz,t,beta1,beta2,n,tau,cov1,cov2)*Sim_deriv1(Nlabel,phi,Di,alpha1,delta,i,zzz,t,beta1,beta2,n,tau,lambda,cov1,cov2))
          }}}}
    KK3[t]=sum(KK4)
  }
  InfKK3=sum(KK3,na.rm=T)
  EndKK3=SusKK1+InfKK3
  K5=rep(0,time)
  for(t in 1:time){
    K6=rep(0,n)
    for(i in 1:n){
      for(zzz in 1:zz){
        if (Nlabel[i]==zzz){
          if(tau[i]==(t+1)){
            K6[i]=-
              as.numeric(cov1[i]*cov2[zzz]*(Sim_Exp(Nlabel,phi,Di,alpha1,delta,lambda,i,zzz,t,beta1,beta2,n,tau,cov1,cov2))^2*Sim_deriv2(Nlabel,phi,Di,alpha1,delta,i
                                                                                                                                         ,zzz,t,beta1,beta2,n,tau,lambda,cov1,cov2))
          }}}}
    K5[t]=sum(K6)
  }
  InfK5=sum(K5,na.rm=T)
  InfBeta1Beta2=EndKK3+InfK5
  ###################################### delta ######################
  S10=rep(0,time)
  for(t in 1:time){
    S9=rep(0,n)
    for(i in 1:n){
      for(zzz in 1:zz){
        if(Nlabel[i]==zzz){
          if (tau[i]>(t+1)|tau[i]== 0){
            S9[i]=-
              as.numeric(Sim_ExpDeriv1(Nlabel,phi,Di,alpha1,delta,lambda,i,zzz,t,beta1
                                       ,beta2,n,tau,cov1,cov2)*exp(phi[zzz]))
          }}}}
    S10[t]=sum(S9)
  }
  SusNew10=sum(S10)
  ###################### infectious period ###################
  I12=rep(0,time)
  for(t in 1:time){
    I11=rep(0,n)
    for(i in 1:n){
      for(zzz in 1:zz){
        if(Nlabel[i]==zzz){
          if(tau[i]==(t+1)){
            I11[i]=as.numeric(Sim_ExpDeriv2(Nlabel,phi,Di,alpha1,delta,lambda,i,zzz,t,beta1,beta2,n,tau,cov1,cov2)*Sim_deriv2(Nlabel,phi,Di,alpha1,delta,i,zzz,t,beta1,beta2,n,tau,lambda,cov1,cov2))-
              as.numeric((Sim_ExpDeriv1(Nlabel,phi,Di,alpha1,delta,lambda1,i,zzz,t,beta1,beta2,n,tau,cov1,cov2))^2*Sim_deriv1(Nlabel,phi,Di,alpha1,delta,i,zzz,t,beta1,beta2,n,tau,lambda,cov1,cov2))
          }}}}
    I12[t]=sum(I11)
  }
  InfeNew12=sum(I12,na.rm=T)
  EndNew12=SusNew10+InfeNew12
  InfDelta=EndNew12
  ######################################### alpha & delta####################
  SO8=rep(0,time)
  for(t in 1:time){
    SO7=rep(0,n)
    for(i in 1:n){
      for(zzz in 1:zz){
        if(Nlabel[i]==zzz){
          if (tau[i]>(t+1)|tau[i]== 0){
            SO7[i]=as.numeric(Sim_ExpDeriv1(Nlabel,phi,Di,alpha1,delta,lambda,i,zzz,
                                            t,beta1,beta2,n,tau,cov1,cov2)*exp(phi[zzz]))
          }}}}
    SO8[t]=sum(SO7)
  }
  SusNewO8=sum(SO8)
  ###################### infectious period ###################
  IO10=rep(0,time)
  for(t in 1:time){
    IO9=rep(0,n)
    for(i in 1:n){
      for(zzz in 1:zz){
        if(Nlabel[i]==zzz){
          if(tau[i]==(t+1)){
            IO9[i]=-as.numeric(Sim_ExpDeriv1(Nlabel,phi,Di,alpha1,delta,lambda,i,zzz,t,beta1,beta2,n,tau,cov1,cov2)*Sim_deriv1(Nlabel,phi,Di,alpha1,delta,i,zzz,t,beta1,beta2,n,tau,lambda,cov1,cov2))
            +as.numeric(Sim_ExpDeriv1(Nlabel,phi,Di,alpha1,delta,lambda1,i,zzz,t,beta1,beta2,n,tau,cov1,cov2)*Sim_Exp(Nlabel,phi,Di,alpha1,delta,lambda1,i,zzz,t,beta1,beta2,n,tau,cov1,cov2)*Sim_deriv2(Nlabel,phi,Di,alpha1,delta,i,zzz,t,beta1,beta2,n,tau,lambda,cov1,cov2))
          }}}}
    IO10[t]=sum(IO9)
  }
  InfeNewO10=sum(IO10,na.rm=T)
  InfAlphaDelta=SusNewO8+InfeNewO10
  ####################################beta1 and alpha ################
  UBy1=rep(0,time)
  for(t in 1:time){
    UBy2=rep(0,n)
    for(i in 1:n){
      for(zzz in 1:zz){
        if (Nlabel[i]==zzz){
          if (tau[i]>(t+1)|tau[i]== 0){
            UBy2[i]=-
              as.numeric(cov1[i]*Sim_Exp(Nlabel,phi,Di,alpha1,delta,lambda,i,zzz,t,beta1,beta2,n,tau,cov1,cov2)*exp(phi[zzz]))
          }}}}
    UBy1[t]=sum(UBy2)
  }
  USusB1=sum(UBy1)
  ###################### infectious period ###################
  UBy3=rep(0,time)
  for(t in 1:time){
    UBy4=rep(0,n)
    for(i in 1:n){
      for(zzz in 1:zz){
        if(Nlabel[i]==zzz){
          if(tau[i]==(t+1)){
            UBy4[i]=as.numeric(cov1[i]*Sim_Exp(Nlabel,phi,Di,alpha1,delta,lambda,i
                                               ,zzz,t,beta1,beta2,n,tau,cov1,cov2)*Sim_deriv1(Nlabel,phi,Di,alpha1,delta,i,zzz,t
                                                                                              ,beta1,beta2,n,tau,lambda,cov1,cov2))-
              as.numeric(cov1[i]*Sim_Exp(Nlabel,phi,Di,alpha1,delta,lambda,i,zzz,t,beta1,beta2,n,tau,cov1,cov2)^2*Sim_deriv2(Nlabel,phi,Di,alpha1,delta,i,zzz,t,beta1
                                                                                                                             ,beta2,n,tau,lambda,cov1,cov2))
          }}}}
    UBy3[t]=sum(UBy4)
  }
  UInfB3=sum(UBy3,na.rm=T)
  InfBeta1Alpha=USusB1+UInfB3
  ####################################beta2 and alpha ################
  W1=rep(0,time)
  for(t in 1:time){
    W2=rep(0,n)
    for(i in 1:n){
      for(zzz in 1:zz){
        if (Nlabel[i]==zzz){
          if (tau[i]>(t+1)|tau[i]== 0){
            W2[i]=-
              as.numeric(cov2[zzz]*Sim_Exp(Nlabel,phi,Di,alpha1,delta,lambda,i,zzz,t
                                           ,beta1,beta2,n,tau,cov1,cov2)*exp(phi[zzz]))
          }}}}
    W1[t]=sum(W2)
  }
  SusW1=sum(W1)
  ###################### infectious period ###################
  W3=rep(0,time)
  for(t in 1:time){
    W4=rep(0,n)
    for(i in 1:n){
      for(zzz in 1:zz){
        if(Nlabel[i]==zzz){
          if(tau[i]==(t+1)){
            W4[i]=as.numeric(cov2[zzz]*Sim_Exp(Nlabel,phi,Di,alpha1,delta,lambda,i
                                               ,zzz,t,beta1,beta2,n,tau,cov1,cov2)*Sim_deriv1(Nlabel,phi,Di,alpha1,delta,i,zzz,t
                                                                                              ,beta1,beta2,n,tau,lambda,cov1,cov2))-
              as.numeric(cov2[zzz]*Sim_Exp(Nlabel,phi,Di,alpha1,delta,lambda,i,zzz,t
                                           ,beta1,beta2,n,tau,cov1,cov2)^2*Sim_deriv2(Nlabel,phi,Di,alpha1,delta,i,zzz,t,beta1,beta2,n,tau,lambda,cov1,cov2))
          }}}}
    W3[t]=sum(W4)
  }
  InfW3=sum(W3,na.rm=T)
  InfBeta2Alpha=SusW1+InfW3
  ###################################beta1 and delta###################
  ################## susceptible period ####################
  PS8=rep(0,time)
  for(t in 1:time){
    PS7=rep(0,n)
    for(i in 1:n){
      for(zzz in 1:zz){
        if(Nlabel[i]==zzz){
          if (tau[i]>(t+1)|tau[i]== 0){
            PS7[i]=as.numeric(cov1[i]*Sim_ExpDeriv1(Nlabel,phi,Di,alpha1,delta,lambda,i,zzz,t,beta1,beta2,n,tau,cov1,cov2)*exp(phi[zzz]))
          }}}}
    PS8[t]=sum(PS7)
  }
  PSusNew8=sum(PS8)
  ###################### infectious period ###################
  PIT10=rep(0,time)
  for(t in 1:time){
    PIT9=rep(0,n)
    for(i in 1:n){
      for(zzz in 1:zz){
        if(Nlabel[i]==zzz){
          if(tau[i]==(t+1)){
            PIT9[i]=-
              as.numeric(cov1[i]*Sim_ExpDeriv1(Nlabel,phi,Di,alpha1,delta,lambda,i,zzz,t,beta1,beta2,n,tau,cov1,cov2)*Sim_deriv1(Nlabel,phi,Di,alpha1,delta,i,zzz,t,beta1,beta2,n,tau,lambda,cov1,cov2))
            +as.numeric(cov1[i]*Sim_ExpDeriv1(Nlabel,phi,Di,alpha1,delta,lambda,i,zzz,t,beta1,beta2,n,tau,cov1,cov2)*Sim_Exp(Nlabel,phi,Di,alpha1,delta,lambda1,i,zzz,t,beta1,beta2,n,tau,cov1,cov2)*Sim_deriv2(Nlabel,phi,Di,alpha1,delta,i,zzz,t,beta1,beta2,n,tau,lambda,cov1,cov2))
          }}}}
    PIT10[t]=sum(PIT9)
  }
  PInfeNew10=sum(PIT10,na.rm=T)
  InfBeta1Delta=PSusNew8+PInfeNew10
  ###################################beta2 and delta###################
  ################## susceptible period ####################
  Z8=rep(0,time)
  for(t in 1:time){
    Z7=rep(0,n)
    for(i in 1:n){
      for(zzz in 1:zz){
        if(Nlabel[i]==zzz){
          if (tau[i]>(t+1)|tau[i]== 0){
            Z7[i]=as.numeric(cov2[zzz]*Sim_ExpDeriv1(Nlabel,phi,Di,alpha1,delta,lambda,i,zzz,t,beta1,beta2,n,tau,cov1,cov2)*exp(phi[zzz]))
          }}}}
    Z8[t]=sum(Z7)
  }
  SusNewZ8=sum(Z8)
  ###################### infectious period ###################
  Z10=rep(0,time)
  for(t in 1:time){
    Z9=rep(0,n)
    for(i in 1:n){
      for(zzz in 1:zz){
        if(Nlabel[i]==zzz){
          if(tau[i]==(t+1)){
            Z9[i]=-
              as.numeric(cov2[zzz]*Sim_ExpDeriv1(Nlabel,phi,Di,alpha1,delta,lambda,i,zzz,t,beta1,beta2,n,tau,cov1,cov2)*Sim_deriv1(Nlabel,phi,Di,alpha1,delta,i,zzz,t,beta1,beta2,n,tau,lambda,cov1,cov2))
            +as.numeric(cov2[zzz]*Sim_ExpDeriv1(Nlabel,phi,Di,alpha1,delta,lambda,i,zzz,t,beta1,beta2,n,tau,cov1,cov2)*Sim_Exp(Nlabel,phi,Di,alpha1,delta,lambda1,i,zzz,t,beta1,beta2,n,tau,cov1,cov2)*Sim_deriv2(Nlabel,phi,Di,alpha1,delta,i,zzz,t,beta1,beta2,n,tau,lambda,cov1,cov2))
          }}}}
    Z10[t]=sum(Z9)
  }
  InfeNewZ10=sum(Z10,na.rm=T)
  InfBeta2Delta=SusNewZ8+InfeNewZ10
  ################################### sigma^2 ######################
  B=(lambda1*D+(1-lambda1)*I)
  InfSigma=-zz/sigma1^2-Sim_Estsigma(phi,D,lambda1,I)
  InfLambda=-tr((D-I)%*%solve(B)%*%(D-I)%*%solve(B))/2
  ################################### sigma^2 & lambda ######################
  InfSigmaLambda=sigma1*Sim_Estlambda(phi,D,I)
  result=list(InfBeta1Beta2=InfBeta1Beta2,InfBeta1Delta=InfBeta1Delta,
              InfBeta2Delta=InfBeta2Delta,InfBeta1Alpha=InfBeta1Alpha,InfBeta2Alpha=InfBeta2Alpha,InfBeta1=InfBeta1,InfBeta2=InfBeta2,InfAlpha=InfAlpha,InfDelta=InfDelta,InfAlphaDelta=InfAlphaDelta,InfSigma=InfSigma,InfLambda=InfLambda,InfSigmaLambda=InfSigmaLambda)
}


for(iter in 1:ITER){
  mu=rep(0,zz)
  Sigma0=solve(sigma0^2*(lambda0*D+(1-lambda0)*I))
  GePHI[iter,]=mvrnorm(1, mu, Sigma0, tol = 1e-6)
  if(iter>s1){
    #####################
    ###########################
    x11=matrix(0,nSample,g)
    for(i in 1:g){
      for(j in 1:8){
        x11[j,i]=runif(1,0,d*1)
      }
      for(j in 9:10){
        x11[j,i]=runif(1,2.9,d*1)
      }
    }
    x1=c(x11)
    x22=matrix(0,nSample,g)
    for(i in 1:g){
      for(j in 1:6){
        x22[j,i]=runif(1,d*1,d*2)
      }
      for(j in 7:8){
        x22[j,i]=runif(1,d*1,d*1+.1)
      }
      for(j in 9:10){
        x22[j,i]=runif(1,d*2-.1,d*2)
      }
    }
    x2=c(x22)
    x33=matrix(0,nSample,g)
    for(i in 1:g){
      for(j in 1:6){
        x33[j,i]=runif(1,d*2,d*3)
      }
      for(j in 7:8){
        x33[j,i]=runif(1,d*2,d*2+.1)
      }
      for(j in 9:10){
        x33[j,i]=runif(1,d*3-.1,d*3)
      }
    }
    x3=c(x33)
    x44=matrix(0,nSample,g)
    for(i in 1:g){
      for(j in 1:6){
        x44[j,i]=runif(1,d*3,d*4)
      }
      for(j in 7:8){
        x44[j,i]=runif(1,d*3,d*3+.1)
      }
      for(j in 9:10){
        x44[j,i]=runif(1,d*4-.1,d*4)
      }
    }
    x4=c(x44)
    x55=matrix(0,nSample,g)
    for(i in 1:g){
      for(j in 1:6){
        x55[j,i]=runif(1,d*4,d*5)
      }
      for(j in 7:8){
        x55[j,i]=runif(1,d*4,d*4+.1)
      }
      for(j in 9:10){
        x55[j,i]=runif(1,d*5-.1,d*5)
      }
    }
    x5=c(x55)
    x66=matrix(0,nSample,g)
    for(i in 1:g){
      for(j in 1:6){
        x66[j,i]=runif(1,d*5,d*6)
      }
      for(j in 7:8){
        x66[j,i]=runif(1,d*5,d*5+.1)
      }
      for(j in 9:10){
        x66[j,i]=runif(1,d*6-.1,d*6)
      }
    }
    x6=c(x66)
    x77=matrix(0,nSample,g)
    for(i in 1:g){
      for(j in 1:6){
        x77[j,i]=runif(1,d*6,d*7)
      }
      for(j in 7:8){
        x77[j,i]=runif(1,d*6,d*6+.1)
      }
      for(j in 9:10){
        x77[j,i]=runif(1,d*7-.1,d*7)
      }
    }
    x7=c(x77)
    x88=matrix(0,nSample,g)
    for(i in 1:g){
      for(j in 1:6){
        x88[j,i]=runif(1,d*7,d*8)
      }
      for(j in 7:8){
        x88[j,i]=runif(1,d*7,d*7+.1)
      }
      for(j in 9:10){
        x88[j,i]=runif(1,d*8-.1,d*8)
      }
    }
    x8=c(x88)
    y11=matrix(0,nSample,g)
    for(i in 1:g){
      for(j in 1:3){
        y11[j,i]=runif(1,d*(i-1),d*i-2.9)
      }
      for(j in 4:7){
        y11[j,i]=runif(1,d*i-.1,d*i)
      }
      for(j in 8:10){
        y11[j,i]=runif(1,d*(i-1),d*i)
      }}
    y1=c(y11)
    y22=matrix(0,nSample,g)
    for(i in 1:g){
      for(j in 1:3){
        y22[j,i]=runif(1,d*(i-1),d*i-2.9)
      }
      for(j in 4:7){
        y22[j,i]=runif(1,d*i-.1,d*i)
      }
      for(j in 8:10){
        y22[j,i]=runif(1,d*(i-1),d*i)
      }}
    y2=c(y22)
    y33=matrix(0,nSample,g)
    for(i in 1:g){
      for(j in 1:3){
        y33[j,i]=runif(1,d*(i-1),d*i-2.9)
      }
      for(j in 4:7){
        y33[j,i]=runif(1,d*i-.1,d*i)
      }
      for(j in 8:10){
        y33[j,i]=runif(1,d*(i-1),d*i)
      }}
    y3=c(y33)
    y44=matrix(0,nSample,g)
    for(i in 1:g){
      for(j in 1:3){
        y44[j,i]=runif(1,d*(i-1),d*i-2.9)
      }
      for(j in 4:7){
        y44[j,i]=runif(1,d*i-.1,d*i)
      }
      for(j in 8:10){
        y44[j,i]=runif(1,d*(i-1),d*i)
      }}
    y4=c(y44)
    y55=matrix(0,nSample,g)
    for(i in 1:g){
      for(j in 1:3){
        y55[j,i]=runif(1,d*(i-1),d*i-2.9)
      }
      for(j in 4:7){
        y55[j,i]=runif(1,d*i-.1,d*i)
      }
      for(j in 8:10){
        y55[j,i]=runif(1,d*(i-1),d*i)
      }}
    y5=c(y55)
    y66=matrix(0,nSample,g)
    for(i in 1:g){
      for(j in 1:3){
        y66[j,i]=runif(1,d*(i-1),d*i-2.9)
      }
      for(j in 4:7){
        y66[j,i]=runif(1,d*i-.1,d*i)
      }
      for(j in 8:10){
        y66[j,i]=runif(1,d*(i-1),d*i)
      }}
    y6=c(y66)
    y77=matrix(0,nSample,g)
    for(i in 1:g){
      for(j in 1:3){
        y77[j,i]=runif(1,d*(i-1),d*i-2.9)
      }
      for(j in 4:7){
        y77[j,i]=runif(1,d*i-.1,d*i)
      }
      for(j in 8:10){
        y77[j,i]=runif(1,d*(i-1),d*i)
      }}
    y7=c(y77)
    y88=matrix(0,nSample,g)
    for(i in 1:g){
      for(j in 1:3){
        y88[j,i]=runif(1,d*(i-1),d*i-2.9)
      }
      for(j in 4:7){
        y88[j,i]=runif(1,d*i-.1,d*i)
      }
      for(j in 8:10){
        y88[j,i]=runif(1,d*(i-1),d*i)
      }}
    y8=c(y88)
    xxx=c(x1,x2,x3,x4,x5,x6,x7,x8)
    yyy=c(y1,y2,y3,y4,y5,y6,y7,y8)
    data=cbind(xxx,yyy)
    x=data[,1]
    y=data[,2]
    Di=matrix(0,n,n)
    for(i in 1:n){
      for(j in 1:n){
        Di[i,j]=sqrt((x[i]-x[j])^2+(y[i]-y[j])^2)
      }}
    Di=Di
    phi=GePHI[iter,]
    Rate=array(0,c(n,time,zz))
    for(t in 1:time){
      for(i in 1:n){
        for(zzz in 1:zz){
          if (Nlabel[i]==zzz){
            if(tau[i]== 0){
              dx1=rep(0,n)
              for(j in 1:n){
                if (Nlabel[j]==zzz){
                  if(tau[j]<=t & (tau[j]+lambda[j])>t & tau[j]!=0){
                    dx1[j]=Di[i,j]^(-delta0)
                  }}}
              dx=sum(dx1)
              Rate[i,t,zzz]=(exp(alpha0+beta10*cov1[i]+beta20*cov2[zzz]
                                 +phi[zzz])*dx)
              P=1-exp(-exp(alpha0+beta10*cov1[i]+beta20*cov2[zzz]+phi[zzz])*dx)
              u=runif(1,0,1)
              if(P>u){
                tau[i]=t+1
              }}}}}}
    tau
    OUT2=list()
    ss=1
    LA=numeric()
    est0=Sim_Estpar(Nlabel,phi,Di,alpha0,delta0,lambda0,sigma0,beta10,beta20,zz,time,n,tau,lambda,I,D,cov1,cov2)
    alpha1=est0$alpha1
    alpha1
    delta=est0$delta
    delta
    lambda1=est0$lambda1
    lambda1
    sigma1=est0$sigma1
    sigma1
    beta1=est0$beta1
    beta1
    beta2=est0$beta2
    beta2
    Uhat=est0$Uhat
    LA[1]=Sim_Loglik(Nlabel,phi,Di,alpha1,delta,lambda,sigma1,beta1,beta2,time,n,zz,tau,lambda1,I,D,cov1,cov2)
    OUT2[[1]]=c(alpha1,beta1,beta2,delta,sigma1,lambda1)
    repeat{
      ss=ss+1
      if(ss>15){
        alpha1=NA
        delta=NA
        sigma1=NA
        lambda1=NA
        out1=list(alpha1=NA,delta=NA,sigma1=NA,lambda1=NA,beta1=NA,beta2=NA,
                  ss=NA)
        break
      }
      if(!is.na(alpha1)&&!is.na(delta)&&!is.na(lambda1)&&!is.na(sigma1)&&!
         is.na(beta1)&&!is.na(beta2)){
        if(sigma1<0 |lambda1>1|lambda1<0){
          alpha1=NA
          delta=NA
          sigma1=NA
          lambda1=NA
          out1=list(alpha1=NA,delta=NA,sigma1=NA,lambda1=NA,beta1=NA,beta2=NA,
                    ss=ss)
          break
        }}
      E0=c(alpha1,beta1,beta2,delta,sigma1,lambda1)
      est=Sim_Estpar(Nlabel,phi,Di,alpha1,delta,lambda1,sigma1,beta1,beta2,zz,time,n,tau,lambda,I,D,cov1,cov2)
      alpha1=est$alpha1
      beta1=est$beta1
      beta2=est$beta2
      delta=est$delta
      lambda1=est$lambda1
      sigma1=est$sigma1
      Uhat=est$Uhat
      if(!is.na(alpha1)&&!is.na(delta)&&!is.na(lambda1)&&!is.na(sigma1)&&!
         is.na(beta1)&&!is.na(beta2)){
        if(sigma1<0 |lambda1>1|lambda1<0){
          alpha1=NA
          delta=NA
          sigma1=NA
          lambda1=NA
          out1=list(alpha1=NA,delta=NA,sigma1=NA,lambda1=NA,beta1=NA,beta2=NA,
                    ss=ss,LA=LA)
          break
        }
        if(sigma1>=0 & lambda1<=1 & lambda1>=0){
          E1=c(alpha1,beta1,beta2,delta,sigma1,lambda1)
          Dist=sqrt(sum((E0-E1)^2))
          LA[ss]=Sim_Loglik(Nlabel,phi,Di,alpha1,delta,lambda,sigma1,beta1,beta2,time,n,zz,tau,lambda1,I,D,cov1,cov2)
          out1=list(Uhat=Uhat,alpha1=alpha1,beta1=beta1,beta2=beta2,delta=delta,sigma1=sigma1,lambda1=lambda1,ss=ss,LA=LA)
          OUT2[[ss]]=c(alpha1,beta1,beta2,delta,sigma1,lambda1)
          if(!Dist>0.5){
            break
          }}}
    }
    out1
    VecEM=c(out1$alpha1,out1$beta1,out1$beta2,out1$delta,out1$sigma1,out1$lambda1)
    Rate1=array(0,c(n,time,zz))
    for(t in 1:time){
      for(i in 1:n){
        for(zzz in 1:zz){
          if (Nlabel[i]==zzz){
            if(tau[i]==(t+1)){
              dx1=rep(0,n)
              for(j in 1:n){
                if (Nlabel[j]==zzz){
                  if(tau[j]<=t & (tau[j]+lambda[j])>t & tau[j]!=0){
                    dx1[j]=Di[i,j]^(-delta)
                  }}}
              dx=sum(dx1)
              Rate1[i,t,zzz]=(exp(out1$alpha1+out1$beta1*cov1[i]
                                  +out1$beta2*cov2[zzz]+phi[zzz])*dx)
            }}}}}
    if(is.na(sigma1)){
      Estpar=c(NA,NA,NA,NA,NA)
      SimPAR[iter-s1,]=NA
      SimRate[,,,iter-s1]=NA
      SimRate1[,,,iter-s1]=NA
      SS[iter]=NA
      ComInf[,,iter-s1]=NA
      Coverage[iter-s1,]=NA
    }
    if(!is.na(sigma1)){
      Estpar=c(alpha1,beta1,beta2,delta,sigma1,lambda1)
      EstInf = Sim_Information(Nlabel,Uhat,Di,VecEM[1],VecEM[2],VecEM[3],VecEM[4],VecEM[5],VecEM[6],time,n,zz,tau,D,cov1,cov2,lambda,I)
      a1=EstInf$InfAlpha
      a2=EstInf$InfBeta1Alpha
      a3=EstInf$InfBeta2Alpha
      a4=EstInf$InfAlphaDelta
      a5=EstInf$InfBeta1
      a6=EstInf$InfBeta1Beta2
      a7=EstInf$InfBeta1Delta
      a8=EstInf$InfBeta2
      a9=EstInf$InfBeta2Delta
      a10=EstInf$InfDelta
      a11=EstInf$InfSigma
      a12=EstInf$InfSigmaLambda
      a13=EstInf$InfLambda
      L1=matrix(c(a1,a2,a3,a4,a2,a5,a6,a7,a3,a6,a8,a9,a4,a7,a9,a10),4,4)
      L2=matrix(c(a11,a12,a12,a13),2,2)
      CInf2=sqrt(diag(-solve(as.matrix(bdiag(L1,L2)))))
      left1=Estpar-CInf2*qnorm(0.975)
      right1=Estpar+CInf2*qnorm(0.975)
      True=c(alpha0,beta10,beta20,delta0,sigma0,lambda0)
      Cover2=rep(0,6)
      for(j in 1:6){
        if(!is.na(left1[1])){
          if(left1[j]<True[j] & right1[j]>True[j]){
            Cover2[j]=1
          }}}
      SimPAR[iter-s1,]=Estpar
      ComInf[,,iter-s1]=-solve(as.matrix(bdiag(L1,L2)))
      SimRate[,,,iter-s1]=Rate
      SimRate1[,,,iter-s1]=Rate1
      Coverage[iter-s1,]=Cover2
      SS[iter-s1]=ss
      MAX[iter-s1]=max(SimRate1[,,,iter-s1])
    }
    ComInf=ComInf
    NewPar=SimPAR
    NewRate1=SimRate1
    NewMAX=MAX
    SS=SS
  }
}
}
