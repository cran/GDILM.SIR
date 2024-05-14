#' Calculating the estimated values for the parameters using log-likelihood function
#'
#' @param Nlabel Label for each sample from the area
#' @param phi Spatial random effects
#' @param Di Euclidean distance between susceptible individual and infectious individual
#' @param alpha1 intercept
#' @param delta Spatial parameter
#' @param lambda1 Spatial dependence
#' @param sigma1 precision of spatial random effects
#' @param beta1 the parameter corresponding to the covariate associated with susceptible individual
#' @param beta2 the parameter corresponding to the covariate associated with area
#' @param zz Number of areas
#' @param time Time
#' @param n Total number of individuals
#' @param tau the set of infectious individuals at time t in the zth area
#' @param lambda a vector containing the length of infectious period
#' @param I identity matrix
#' @param D Neighborhood structure
#' @param cov1 Individual level covariates
#' @param cov2 Area level covariates
#'
#' @return a list of the solutions for the estimations of the parameters
#' @export
#'
#' @importFrom stats optim
#' @importFrom mvtnorm dmvnorm
#' @importFrom psych tr
#'
#'
#' @examples Sim_Estpar(rep(1:4,each=5),runif(4,min = 0, max = 1),
#' matrix(runif(400,min=4,max=20),nrow=20,byrow = TRUE),0.4,3,0.2,0.5,1,1,4,10,
#' 20,sample(c(0,1),replace = TRUE, size = 20),rep(3,20),diag(4),
#' matrix(c(0,-1,0,-1,-1,0,-1,-1,0,-1,0,-1,-1,-1,-1,0),nrow=4,byrow=TRUE),
#' runif(20, 0, 1),runif(4, 0, 1))
#'
#'
#'
Sim_Estpar=function(Nlabel,phi,Di,alpha1,delta,lambda1,sigma1,beta1,beta2,zz,
                    time,n,tau,lambda,I,D,cov1,cov2){
####################################

  Loglik1=function(par){

    FL1=rep(0,time)
    for(t in 1:time){
      FL2=rep(0,n)
      for(i in 1:n){
        for(zzz in 1:zz){
          if(Nlabel[i]==zzz){
            if(tau[i]>(t+1)|tau[i]==0){
              dx1=rep(0,n)
              for(j in 1:n){
                if (Nlabel[j]==zzz){
                  if(tau[j]<=t & (tau[j]+lambda[j])>t & tau[j]!=0){
                    dx1[j]=Di[i,j]^(-delta)
                  }}}
              dx1=replace(dx1,is.infinite(dx1),0)
              dx=sum(dx1)
              prob1=1-exp(-exp(alpha1+beta1*cov1[i]+beta2*cov2[zzz]+phi[zzz])*dx)
              FL2[i]=log(1-prob1)
              FL2[i]=replace(FL2[i],is.infinite(FL2[i]),0)
            }
            if(tau[i]==(t+1)){
              dx4=rep(0,n)
              for(j in 1:n){
                if (Nlabel[j]==zzz){
                  if(tau[j]<=t & (tau[j]+lambda[j])>t & tau[j]!=0){
                    dx4[j]=Di[i,j]^(-delta)
                  }}}
              dx4=replace(dx4,is.infinite(dx4),0)
              dx5=sum(dx4)
              prob1=1-exp(-exp(alpha1+beta1*cov1[i]+beta2*cov2[zzz]+phi[zzz])*dx5)
              FL2[i]=log(prob1)
              FL2[i]=replace(FL2[i],is.infinite(FL2[i]),0)
            }}}}
      FL1[t]=sum(FL2,na.rm=T)
    }
    Sigma1 = solve(sigma1^2*(lambda1*D+(1-lambda1)*I))
    LLL=dmvnorm(phi, rep(0,zz),
                solve(Sigma1^2*(as.numeric(lambda1)*D+(1-as.numeric(lambda1))*I)),log = TRUE)
    loglik1=-((sum(FL1,na.rm=T))+LLL)
  }
  init =as.integer(c(phi))

  ml <- optim( init ,Loglik1)

  EstU=ml$par

  ##############################################################################################

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


  ######################################## alpha #######################################
  A1=rep(0,time)
  for(t in 1:time){
    A2=rep(0,n)
    for(i in 1:n){
      for(zzz in 1:zz){
        if (Nlabel[i]==zzz){
          if (tau[i]>(t+1)|tau[i]== 0){
            A2[i]=-as.numeric(Sim_Exp(Nlabel,phi,Di,alpha1,delta,lambda,i,zzz,t,beta1,beta2,n,tau,cov1,cov2)*exp(EstU[zzz]))
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
            A4[i]=as.numeric(Sim_Exp(Nlabel,phi,Di,alpha1,delta,lambda,i,zzz,t,beta1,beta2,n,tau,cov1,cov2)*Sim_deriv1(Nlabel,EstU,Di,alpha1,delta,i,zzz,t,beta1,beta2,n,tau,lambda,cov1,cov2))
          }}}}
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
            A6[i]=-as.numeric((Sim_Exp(Nlabel,phi,Di,alpha1,delta,lambda,i,zzz,t,beta1,beta2,n,tau,cov1,cov2))^2*Sim_deriv2(Nlabel,EstU,Di,alpha1,delta,i,zzz,t,beta1,beta2,n,tau,lambda,cov1,cov2))
          }}}}
    A5[t]=sum(A6)
  }
  InfA5=sum(A5,na.rm=T)
  EndA5=EndA3+InfA5
  EstAlpha=alpha1-EndA3/EndA5 #solution for alpha

  #################################beta1##################################
  B1=rep(0,time)
  for(t in 1:time){
    B2=rep(0,n)
    for(i in 1:n){
      for(zzz in 1:zz){
        if (Nlabel[i]==zzz){
          if (tau[i]>(t+1)|tau[i]== 0){
            B2[i]=-
              as.numeric(cov1[i]*Sim_Exp(Nlabel,phi,Di,EstAlpha,delta,lambda,i,zzz,t,beta1,beta2,n,tau,cov1,cov2)*exp(EstU[zzz]))
          }}}}
    B1[t]=sum(B2)
  }
  SusB1=sum(B1)
  ###################### infectious period ###################
  B3=rep(0,time)
  for(t in 1:time){
    B4=rep(0,n)
    for(i in 1:n){
      for(zzz in 1:zz){
        if(Nlabel[i]==zzz){
          if(tau[i]==(t+1)){
            B4[i]=as.numeric(cov1[i]*Sim_Exp(Nlabel,phi,Di,EstAlpha,delta,lambda,i,zzz,t,beta1,beta2,n,tau,cov1,cov2)*Sim_deriv1(Nlabel,EstU,Di,alpha1,delta,i,zzz,t,beta1,beta2,n,tau,lambda,cov1,cov2))
          }}}}
    B3[t]=sum(B4)
  }
  InfB3=sum(B3,na.rm=T)
  EndB3=SusB1+InfB3
  ######################### Second Der...########################
  ##################### infectious period ###################
  BB1=rep(0,time)
  for(t in 1:time){
    BB2=rep(0,n)
    for(i in 1:n){
      for(zzz in 1:zz){
        if (Nlabel[i]==zzz){
          if (tau[i]>(t+1)|tau[i]== 0){
            BB2[i]=-
              as.numeric(cov1[i]^2*Sim_Exp(Nlabel,phi,Di,EstAlpha,delta,lambda,i,zzz,t,beta1,beta2,n,tau,cov1,cov2)*exp(EstU[zzz]))
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
            BB4[i]=as.numeric(cov1[i]^2*Sim_Exp(Nlabel,phi,Di,EstAlpha,delta,lambda,i,zzz,t,beta1,beta2,n,tau,cov1,cov2)*Sim_deriv1(Nlabel,EstU,Di,alpha1,delta,i,zzz,t,beta1,beta2,n,tau,lambda,cov1,cov2))
          }
        }}}
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
              as.numeric(cov1[i]^2*(Sim_Exp(Nlabel,phi,Di,EstAlpha,delta,lambda,i,zzz,t,beta1,beta2,n,tau,cov1,cov2))^2*Sim_deriv2(Nlabel,EstU,Di,alpha1,delta,i,zzz,t,beta1,beta2,n,tau,lambda,cov1,cov2))
          }}}}
    B5[t]=sum(B6)
  }
  InfB5=sum(B5,na.rm=T)
  EndB5=EndBB3+InfB5
  EstBeta1=beta1-EndB3/EndB5 # solution for beta1

  #################################beta2##################################
  Y1=rep(0,time)
  for(t in 1:time){
    Y2=rep(0,n)
    for(i in 1:n){
      for(zzz in 1:zz){
        if (Nlabel[i]==zzz){
          if (tau[i]>(t+1)|tau[i]== 0){
            Y2[i]=-
              as.numeric(cov2[zzz]*Sim_Exp(Nlabel,phi,Di,EstAlpha,delta,lambda,i,zzz,t,EstBeta1,beta2,n,tau,cov1,cov2)*exp(EstU[zzz]))
          }}}}
    Y1[t]=sum(Y2)
  }
  SusY1=sum(Y1)
  ###################### infectious period ###################
  Y3=rep(0,time)
  for(t in 1:time){
    Y4=rep(0,n)
    for(i in 1:n){
      for(zzz in 1:zz){
        if(Nlabel[i]==zzz){
          if(tau[i]==(t+1)){
            Y4[i]=as.numeric(cov2[zzz]*Sim_Exp(Nlabel,phi,Di,EstAlpha,delta,lambda,i,zzz,t,EstBeta1,beta2,n,tau,cov1,cov2)*Sim_deriv2(Nlabel,EstU,Di,alpha1,delta,i,zzz,t,beta1,beta2,n,tau,lambda,cov1,cov2))
          }}}}
    Y3[t]=sum(Y4)
  }
  InfY3=sum(Y3,na.rm=T)
  EndY3=SusY1+InfY3
  ######################### Second Der...########################
  ##################### infectious period ###################
  YY1=rep(0,time)
  for(t in 1:time){
    YY2=rep(0,n)
    for(i in 1:n){
      for(zzz in 1:zz){
        if (Nlabel[i]==zzz){
          if (tau[i]>(t+1)|tau[i]== 0){
            YY2[i]=-
              as.numeric(cov2[zzz]^2*Sim_Exp(Nlabel,phi,Di,EstAlpha,delta,lambda,i,zzz,t,EstBeta1,beta2,n,tau,cov1,cov2)*exp(EstU[zzz]))
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
            YY4[i]=as.numeric(cov2[zzz]^2*Sim_Exp(Nlabel,phi,Di,EstAlpha,delta,lambda,i,zzz,t,EstBeta1,beta2,n,tau,cov1,cov2)*Sim_deriv1(Nlabel,EstU,Di,alpha1,delta,i,zzz,t,beta1,beta2,n,tau,lambda,cov1,cov2))
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
              as.numeric(cov2[zzz]^2*(Sim_Exp(Nlabel,phi,Di,EstAlpha,delta,lambda,i,zzz,t,EstBeta1,beta2,n,tau,cov1,cov2))^2*Sim_deriv2(Nlabel,EstU,Di,alpha1,delta,i,zzz,t,beta1,beta2,n,tau,lambda,cov1,cov2))
          }}}}
    Y5[t]=sum(Y6)
  }
  InfY5=sum(Y5,na.rm=T)
  EndY5=EndYY3+InfY5
  EstBeta2=beta2-EndY3/EndY5 # solution for beta2

  ################################# Delta #############################
  ################## susceptible period ####################
  S8=rep(0,time)
  for(t in 1:time){
    S7=rep(0,n)
    for(i in 1:n){
      for(zzz in 1:zz){
        if(Nlabel[i]==zzz){
          if (tau[i]>(t+1)|tau[i]== 0){
            S7[i]=as.numeric(Sim_ExpDeriv1(Nlabel,phi,Di,EstAlpha,delta,lambda,i,zzz,t,EstBeta1,EstBeta2,n,tau,cov1,cov2)*exp(EstU[zzz]))
          }}}}
    S8[t]=sum(S7)
  }
  SusNew8=sum(S8)
  ###################### infectious period ###################
  I10=rep(0,time)
  for(t in 1:time){
    I9=rep(0,n)
    for(i in 1:n){
      for(zzz in 1:zz){
        if(Nlabel[i]==zzz){
          if(tau[i]==(t+1)){
            I9[i]=-
              as.numeric(Sim_ExpDeriv1(Nlabel,phi,Di,EstAlpha,delta,lambda,i,zzz,t,EstBeta1,EstBeta2,n,tau,cov1,cov2)*Sim_deriv1(Nlabel,EstU,Di,alpha1,delta,i,zzz,t,beta1,beta2,n,tau,lambda,cov1,cov2))
          }}}}
    I10[t]=sum(I9)
  }
  InfeNew10=sum(I10,na.rm=T)
  EndNew10=SusNew8+InfeNew10
  ######################### Second Der...########################
  S10=rep(0,time)
  for(t in 1:time){
    S9=rep(0,n)
    for(i in 1:n){
      for(zzz in 1:zz){
        if(Nlabel[i]==zzz){
          if (tau[i]>(t+1)|tau[i]== 0){
            S9[i]=-as.numeric(Sim_ExpDeriv2(Nlabel,phi,Di,EstAlpha,delta,lambda,i,zzz,t,EstBeta1,EstBeta2,n,tau,cov1,cov2)*exp(EstU[zzz]))
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
            I11[i]=as.numeric(Sim_ExpDeriv2(Nlabel,phi,Di,EstAlpha,delta,lambda,i,zzz,t,EstBeta1,EstBeta2,n,tau,cov1,cov2)*Sim_deriv1(Nlabel,EstU,Di,alpha1,delta,i,zzz,t,beta1,beta2,n,tau,lambda,cov1,cov2))-as.numeric((Sim_ExpDeriv1(Nlabel,phi,Di,EstAlpha,delta,lambda,i,zzz,t,EstBeta1,EstBeta2,n,tau,cov1,cov2))^2*Sim_deriv2(Nlabel,EstU,Di,alpha1,delta,i,zzz,t,beta1,beta2,n,tau,lambda,cov1,cov2))
          }}}}
    I12[t]=sum(I11)
  }
  InfeNew12=sum(I12,na.rm=T)
  EndNew12=SusNew10+InfeNew12
  Estdelta=delta-EndNew10/EndNew12 # solution for delta

  ##################################### Sigma_U #####################

  B=(lambda1*D+(1-lambda1)*I)
  firstS=zz/sigma1-sigma1*Sim_Estsigma(EstU,D,lambda1,I)
  firstL=tr((D-I)%*%solve(B))/2+sigma1^2*Sim_Estlambda(EstU,D,I)/2
  firstLS=sigma1*Sim_Estlambda(EstU,D,I)
  secondS=-zz/sigma1^2-Sim_Estsigma(EstU,D,lambda1,I)
  secondL=-tr((D-I)%*%solve(B)%*%(D-I)%*%solve(B))/2
  Inf1=matrix(c(secondS,firstLS,firstLS,secondL),2,2)
  Score=c(firstS,firstL)
  c(sigma1,lambda1)-solve(Inf1)%*%Score
  Aest=c(sigma1,lambda1)-solve(Inf1)%*%Score
  HatSigmmaU=Aest[1] # solution for tau
  EstGammau=Aest[2] # solution for gamma
  ##############################################
  result=list(beta1=EstBeta1,beta2=EstBeta2,Uhat=EstU,alpha1=EstAlpha,
              delta=Estdelta,sigma1=HatSigmmaU,lambda1=EstGammau)
  result
}
