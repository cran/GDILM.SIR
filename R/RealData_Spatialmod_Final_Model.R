#' This function is used to estimate model parameters
#'
#' @param ITER Number of iterations
#' @param zz Number of Regions
#' @param lambda0 Spatial dependence
#' @param sigma0 precision
#' @param Di Euclidean distance between susceptible individual and infectious individual
#' @param D Neighborhood structure
#' @param n total number of individuals
#' @param time time
#' @param tau tau
#' @param lambda lambda ###
#' @param alpha0 intercept
#' @param q1 Number of variables corresponding to individual level data
#' @param q2 Number of variables corresponding to area level data
#' @param cov1 Individual level covariates
#' @param cov2 Area level covariates
#' @param phi Spatial random effects
#' @param delta0 Spatial parameter
#' @param Nlabel Label for each sample from the area
#' @param npar number of parameters
#' @param I Identity matrix
#'
#'
#' @return Numerical values for estimates
#' @export
#'
#' @importFrom MASS mvrnorm
#' @importFrom numDeriv hessian
#' @importFrom stats optim
#' @importFrom mvtnorm dmvnorm
#' @importFrom psych tr
#'
#' @examples Realdata_Finalmodel(2,4,0.2,0.5,
#' matrix(runif(400,min = 4,max = 20),nrow=20, byrow = TRUE),
#' matrix(c(0,-1,0,-1,-1,0,-1,-1,0,-1,0,-1,-1,-1,-1,0),nrow=4,byrow=TRUE),20,10,
#' sample(c(0,1),replace = TRUE, size = 20),rep(3,20),0.4,6,5,
#' matrix(runif(120, 0, 1),nrow=20,byrow=TRUE),
#' matrix(runif(20, 0, 1),nrow=4,byrow=TRUE),runif(4,min = 0, max = 1),2,
#' rep(1:4,each=5),15,diag(4))
#'
Realdata_Finalmodel <- function(ITER,zz,lambda0,sigma0,Di, D,n,time,tau,lambda,alpha0,q1,q2,cov1,cov2,phi,delta0,Nlabel,npar,I){
  beta10=rep(0,q1)
  beta20=rep(0,q2)
  VecEM = matrix(0,ITER,npar)
  EstLoglikObs = c()
  Finalphi = list()
  phi1 = matrix(0,ITER,zz)
  RealLog = c()
  EstUhat = c()
  DiffValue = c()

  Real_Exp=function(Nlabel,phi,Di,alpha1,delta,lambda,i,zzz,t,beta1,beta2,n,tau,cov1,cov2){
    SumH1=rep(0,n)
    for(j in 1:n){
      if (Nlabel[j]==zzz){
        if(tau[j]<=t & (tau[j]+lambda[j])>t & tau[j]!=0){
          SumH1[j]=Di[i,j]^(-delta)
        }
      }
    }
    SumH1=replace(SumH1,is.infinite(SumH1),0)
    SumH2=sum(SumH1)
    SumH3=exp(alpha1+cov1[i,]%*%beta1+cov2[zzz,]%*%beta2)*SumH2
    return(SumH3)
  }

  Real_ExpDeriv1=function(Nlabel,phi,Di,alpha1,delta,lambda,i,zzz,t,beta1,beta2,n,tau,cov1,cov2){
    SumH4=rep(0,n)
    for(j in 1:n){
      if (Nlabel[j]==zzz){
        if(tau[j]<=t & (tau[j]+lambda[j])>t & tau[j]!=0){
          SumH4[j]=Di[i,j]^(-delta)*(log(Di[i,j]))
        }
      }
    }
    SumH4=replace(SumH4,is.infinite(SumH4),0)
    SumH5=sum(SumH4)
    SumH6=exp(alpha1+cov1[i,]%*%beta1+cov2[zzz,]%*%beta2)*SumH5
    return(SumH6)
  }

  Real_ExpDeriv2=function(Nlabel,phi,Di,alpha1,delta,lambda,i,zzz,t,beta1,beta2,n,tau,cov1,cov2){
    SumH7=rep(0,n)
    for(j in 1:n){
      if (Nlabel[j]==zzz){
        if(tau[j]<=t & (tau[j]+lambda[j])>t & tau[j]!=0){
          SumH7[j]=Di[i,j]^(-delta)*(log(Di[i,j]))^2
        }
      }
    }
    SumH7=replace(SumH7,is.infinite(SumH7),0)
    SumH8=sum(SumH7)
    SumH9=exp(alpha1+cov1[i,]%*%beta1+cov2[zzz,]%*%beta2)*SumH8
    return(SumH9)
  }

  Real_prob=function(Nlabel,phi,Di,alpha1,delta,i,zzz,t,beta1,beta2,n,tau,lambda,cov1,cov2){
    dx1=rep(0,n)
    for(j in 1:n){
      if (Nlabel[j]==zzz){
        if(tau[j]<=t & (tau[j]+lambda[j])>t & tau[j]!=0){
          dx1[j]=Di[i,j]^(-delta)
        }
      }
    }
    dx1=replace(dx1,is.infinite(dx1),0)
    dx=sum(dx1)
    prob1=1-exp(-exp(alpha1+cov1[i,]%*%beta1+cov2[zzz,]
                     %*%beta2+phi[zzz])*dx)
    return(prob1)
  }

  Real_deriv1=function(Nlabel,phi,Di,alpha1,delta,i,zzz,t,beta1,beta2,n,tau,lambda,cov1,cov2){
    (1- Real_prob(Nlabel,phi,Di,alpha1,delta,i,zzz,t,beta1,beta2,n,tau,lambda,cov1,cov2))/
      Real_prob(Nlabel,phi,Di,alpha1,delta,i,zzz,t,beta1,beta2,n,tau,lambda,cov1,cov2)*exp(phi[zzz])
  }

  Real_deriv2=function(Nlabel,phi,Di,alpha1,delta,i,zzz,t,beta1,beta2,n,tau,lambda,cov1,cov2){
    (1- Real_prob(Nlabel,phi,Di,alpha1,delta,i,zzz,t,beta1,beta2,n,tau,lambda,cov1,cov2))/
      Real_prob(Nlabel,phi,Di,alpha1,delta,i,zzz,t,beta1,beta2,n,tau,lambda,cov1,cov2)^2*exp(2*phi[zzz
      ])
  }

  Real_Estlambda=function(phi,D,I){
    as.numeric(t(phi)%*%(I-D)%*%phi) #A(lambda)
  }

  Real_Estsigma=function(phi,D,lambda1,I){
    as.numeric(t(phi)%*%(lambda1*D+(1-lambda1)*I)%*%phi) #A(tau)
  }

  RealData_Estpar=function(Nlabel,phi,Di,alpha1,delta,lambda1,sigma1,beta1,beta2,zz,time,n,tau,lambda,D,I,cov1,cov2,q1,q2){
    ####################################


    Loglik1=function(par){
      #for (zzz in 1:zz) {
      #phi[zzz] = par[zzz]
      #}
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
                    }
                  }
                }
                dx1=replace(dx1,is.infinite(dx1),0)
                dx=sum(dx1)
                prob1=1-exp(-exp(alpha1+cov1[i,]%*%beta1+cov2[zzz,]
                                 %*%beta2+phi[zzz])*dx)
                FL2[i]=log(1-prob1)
                FL2[i]=replace(FL2[i],is.infinite(FL2[i]),0)
              }
              if(tau[i]==(t+1)){
                dx4=rep(0,n)
                for(j in 1:n){
                  if (Nlabel[j]==zzz){
                    if(tau[j]<=t & (tau[j]+lambda[j])>t & tau[j]!=0){
                      dx4[j]=Di[i,j]^(-delta)
                    }
                  }
                }
                dx4=replace(dx4,is.infinite(dx4),0)
                dx5=sum(dx4)
                prob1=1-exp(-exp(alpha1+cov1[i,]%*%beta1+cov2[zzz,]
                                 %*%beta2+phi[zzz])*dx5)
                FL2[i]=log(prob1)
                FL2[i]=replace(FL2[i],is.infinite(FL2[i]),0)
              }
            }
          }
        }
        FL1[t]=sum(FL2,na.rm=T)
      }
      Sigma1 = solve(sigma1^2*(lambda1*D+(1-lambda1)*I))
      LLL=dmvnorm(phi, rep(0,zz),
                  solve(Sigma1^2*(as.numeric(lambda1)*D+(1-as.numeric(lambda1))*I)),
                  log = TRUE)
      loglik1=-((sum(FL1,na.rm=T))+LLL)
    }
    init =as.integer(c(phi))
    #init = c(phi)
    ml <- optim( init ,Loglik1)

    EstU=ml$par

    ######################################## alpha #######################################
    A1=rep(0,time)
    for(t in 1:time){
      A2=rep(0,n)
      for(i in 1:n){
        for(zzz in 1:zz){
          if (Nlabel[i]==zzz){
            if (tau[i]>(t+1)|tau[i]== 0){
              A2[i]=-as.numeric(Real_Exp(Nlabel,phi,Di,alpha1,delta,lambda,i,zzz,t,beta1,beta2,n,tau,cov1,cov2)*exp(EstU[zzz]))
            }
          }
        }
      }
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
              A4[i]=as.numeric(Real_Exp(Nlabel,phi,Di,alpha1,delta,lambda,i,zzz,t,beta1,beta2,n,tau,cov1,cov2)*Real_deriv1(Nlabel,EstU,Di,alpha1,delta,i,zzz,t,beta1,beta2,n,tau,lambda,cov1,cov2))
            }
          }
        }
      }
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
              A6[i]=-as.numeric((Real_Exp(Nlabel,phi,Di,alpha1,delta,lambda,i,zzz,t,beta1,beta2,n,tau,cov1,cov2))^2*Real_deriv2(Nlabel,EstU,Di,alpha1,delta,i,zzz,t,beta1,beta2,n,tau,lambda,cov1,cov2))
            }
          }
        }
      }
      A5[t]=sum(A6)
    }
    InfA5=sum(A5,na.rm=T)
    EndA5=EndA3+InfA5
    EstAlpha=alpha1-EndA3/EndA5
    #################################beta1##############################
    ####
    B1=array(0,c(q1,1,time))
    for(t in 1:time){
      B2=array(0,c(q1,1,n))
      for(i in 1:n){
        for(zzz in 1:zz){
          if (Nlabel[i]==zzz){
            if (tau[i]>(t+1)|tau[i]== 0){
              B2[,,i]=-cov1[i,]*as.numeric(Real_Exp(Nlabel,phi,Di,EstAlpha,delta,lambda,i,zzz,t,beta1,beta2,n,tau,cov1,cov2)*exp(EstU[zzz]))
            }
          }
        }
      }
      B1[,,t]=apply(B2,c(1,2),sum)
    }
    SusB1=apply(B1,c(1,2),sum)
    ###################### infectious period ###################
    B3=array(0,c(q1,1,time))
    for(t in 1:time){
      B4=array(0,c(q1,1,n))
      for(i in 1:n){
        for(zzz in 1:zz){
          if(Nlabel[i]==zzz){
            if(tau[i]==(t+1)){
              B4[,,i]=cov1[i,]*as.numeric(Real_Exp(Nlabel,phi,Di,EstAlpha,delta,lambda,i,zzz,t,beta1,beta2,n,tau,cov1,cov2)*Real_deriv1(Nlabel,EstU,Di,alpha1,delta,i,zzz,t,beta1,beta2,n,tau,lambda,cov1,cov2))
            }
          }
        }
      }
      B3[,,t]=apply(B4,c(1,2),sum,na.rm=T)
    }
    InfB3=apply(B3,c(1,2),sum,na.rm=T)
    EndB3=SusB1+InfB3
    ######################### Second Der...########################
    ##################### infectious period ###################
    BB1=array(0,c(q1,q1,time))
    for(t in 1:time){
      BB2=array(0,c(q1,q1,n))
      for(i in 1:n){
        for(zzz in 1:zz){
          if (Nlabel[i]==zzz){
            if (tau[i]>(t+1)|tau[i]== 0){
              BB2[,,i]=-cov1[i,]%*%t(cov1[i,])*as.numeric(Real_Exp(Nlabel,phi,Di,EstAlpha,delta,lambda,i,zzz,t,beta1,beta2,n,tau,cov1,cov2)*exp(EstU[zzz]))
            }
          }
        }
      }
      BB1[,,t]=apply(BB2,c(1,2),sum,na.rm=T)
    }
    SusBB1=apply(BB1,c(1,2),sum,na.rm=T)
    ###################### infectious period ###################
    BB3=array(0,c(q1,q1,time))
    for(t in 1:time){
      BB4=array(0,c(q1,q1,n))
      for(i in 1:n){
        for(zzz in 1:zz){
          if(Nlabel[i]==zzz){
            if(tau[i]==(t+1)){
              BB4[,,i]=cov1[i,]%*%t(cov1[i,])*as.numeric(Real_Exp(Nlabel,phi,Di,EstAlpha,delta,lambda,i,zzz,t,beta1,beta2,n,tau,cov1,cov2)*Real_deriv1(Nlabel,EstU,Di,alpha1,delta,i,zzz,t,beta1,beta2,n,tau,lambda,cov1,cov2))
            }
          }
        }
      }
      BB3[,,t]=apply(BB4,c(1,2),sum,na.rm=T)
    }
    InfBB3=apply(BB3,c(1,2),sum,na.rm=T)
    EndBB3=SusBB1+InfBB3
    B5=array(0,c(q1,q1,time))
    for(t in 1:time){
      B6=array(0,c(q1,q1,n))
      for(i in 1:n){
        for(zzz in 1:zz){
          if (Nlabel[i]==zzz){
            if(tau[i]==(t+1)){
              B6[,,i]=-cov1[i,]%*%t(cov1[i,])*as.numeric((Real_Exp(Nlabel,phi,Di,EstAlpha,delta,lambda,i,zzz,t,beta1,beta2,n,tau,cov1,cov2))^2*Real_deriv2(Nlabel,EstU,Di,alpha1,delta,i,zzz,t,beta1,beta2,n,tau,lambda,cov1,cov2))
            }
          }
        }
      }
      B5[,,t]=apply(B6,c(1,2),sum,na.rm=T)
    }
    InfB5=apply(B5,c(1,2),sum,na.rm=T)
    EndB5=EndBB3+InfB5
    EstBeta1=beta1-solve(EndB5)%*%EndB3
    #################################beta2##############################
    ####
    Y1=array(0,c(q2,1,time))
    for(t in 1:time){
      Y2=array(0,c(q2,1,n))
      for(i in 1:n){
        for(zzz in 1:zz){
          if (Nlabel[i]==zzz){
            if (tau[i]>(t+1)|tau[i]== 0){
              Y2[,,i]=-
                cov2[zzz,]*as.numeric(Real_Exp(Nlabel,phi,Di,EstAlpha,delta,lambda,i,zzz,t,EstBeta1,beta2,n,tau,cov1,cov2)*exp(EstU[zzz]))
            }
          }
        }
      }
      Y1[,,t]=apply(Y2,c(1,2),sum,na.rm=T)
    }
    SusY1=apply(Y1,c(1,2),sum,na.rm=T)
    ###################### infectious period ###################
    Y3=array(0,c(q2,1,time))
    for(t in 1:time){
      Y4=array(0,c(q2,1,n))
      for(i in 1:n){
        for(zzz in 1:zz){
          if(Nlabel[i]==zzz){
            if(tau[i]==(t+1)){
              Y4[,,i]=cov2[zzz,]*as.numeric(Real_Exp(Nlabel,phi,Di,EstAlpha,delta,lambda,i,zzz,t,EstBeta1,beta2,n,tau,cov1,cov2)*Real_deriv1(Nlabel,EstU,Di,alpha1,delta,i,zzz,t,beta1,beta2,n,tau,lambda,cov1,cov2))
            }
          }
        }
      }
      Y3[,,t]=apply(Y4,c(1,2),sum,na.rm=T)
    }
    InfY3=apply(Y3,c(1,2),sum,na.rm=T)
    EndY3=SusY1+InfY3
    ######################### Second Der...########################
    ##################### infectious period ###################
    YY1=array(0,c(q2,q2,time))
    for(t in 1:time){
      YY2=array(0,c(q2,q2,n))
      for(i in 1:n){
        for(zzz in 1:zz){
          if (Nlabel[i]==zzz){
            if (tau[i]>(t+1)|tau[i]== 0){
              YY2[,,i]=-cov2[zzz,]%*%t(cov2[zzz,])*as.numeric(Real_Exp(Nlabel,phi,Di,EstAlpha,delta,lambda,i,zzz,t,EstBeta1,beta2,n,tau,cov1,cov2)*exp(EstU[zzz]))
            }
          }
        }
      }
      YY1[,,t]=apply(YY2,c(1,2),sum,na.rm=T)
    }
    SusYY1=apply(YY1,c(1,2),sum,na.rm=T)
    ###################### infectious period ###################
    YY3=array(0,c(q2,q2,time))
    for(t in 1:time){
      YY4=array(0,c(q2,q2,n))
      for(i in 1:n){
        for(zzz in 1:zz){
          if(Nlabel[i]==zzz){
            if(tau[i]==(t+1)){
              YY4[,,i]=cov2[zzz,]%*%t(cov2[zzz,])*as.numeric(Real_Exp(Nlabel,phi,Di,EstAlpha,delta,lambda,i,zzz,t,EstBeta1,beta2,n,tau,cov1,cov2)*Real_deriv1(Nlabel,EstU,Di,alpha1,delta,i,zzz,t,beta1,beta2,n,tau,lambda,cov1,cov2))
            }
          }
        }
      }
      YY3[,,t]=apply(YY4,c(1,2),sum,na.rm=T)
    }
    InfYY3=apply(YY3,c(1,2),sum,na.rm=T)
    EndYY3=SusYY1+InfYY3
    Y5=array(0,c(q2,q2,time))
    for(t in 1:time){
      Y6=array(0,c(q2,q2,n))
      for(i in 1:n){
        for(zzz in 1:zz){
          if (Nlabel[i]==zzz){
            if(tau[i]==(t+1)){
              Y6[,,i]=-cov2[zzz,]%*%t(cov2[zzz,])*as.numeric((Real_Exp(Nlabel,phi,Di,EstAlpha,delta,lambda,i,zzz,t,EstBeta1,beta2,n,tau,cov1,cov2))^2*Real_deriv2(Nlabel,EstU,Di,alpha1,delta,i,zzz,t,beta1,beta2,n,tau,lambda,cov1,cov2))
            }
          }
        }
      }
      Y5[,,t]=apply(Y6,c(1,2),sum,na.rm=T)
    }
    InfY5=apply(Y5,c(1,2),sum,na.rm=T)
    EndY5=EndYY3+InfY5
    epsilon <- 1e-3
    if (det(EndY5) < epsilon) {
      EndY5 <- EndY5 + diag(epsilon, nrow(EndY5))
    }
    EstBeta2=beta2-solve(EndY5)%*%EndY3
    ################################# Delta
    #############################
    ################## susceptible period ####################
    S8=rep(0,time)
    for(t in 1:time){
      S7=rep(0,n)
      for(i in 1:n){
        for(zzz in 1:zz){
          if(Nlabel[i]==zzz){
            if (tau[i]>(t+1)|tau[i]== 0){
              S7[i]=as.numeric(Real_ExpDeriv1(Nlabel,phi,Di,EstAlpha,delta,lambda,i,zzz,t,EstBeta1,EstBeta2,n,tau,cov1,cov2)*exp(EstU[zzz]))
            }
          }
        }
      }
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
              I9[i]=-as.numeric(Real_ExpDeriv1(Nlabel,phi,Di,EstAlpha,delta,lambda,i,zzz,t,EstBeta1,EstBeta2,n,tau,cov1,cov2)*Real_deriv1(Nlabel,EstU,Di,alpha1,delta,i,zzz,t,beta1,beta2,n,tau,lambda,cov1,cov2))
            }
          }
        }
      }
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
              S9[i]=-as.numeric(Real_ExpDeriv1(Nlabel,phi,Di,EstAlpha,delta,lambda,i,zzz,t,EstBeta1,EstBeta2,n,tau,cov1,cov2)*exp(EstU[zzz]))
            }
          }
        }
      }
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
              I11[i]=as.numeric(Real_ExpDeriv2(Nlabel,phi,Di,EstAlpha,delta,lambda,i,zzz,t,EstBeta1,EstBeta2,n,tau,cov1,cov2)*Real_deriv1(Nlabel,EstU,Di,alpha1,delta,i,zzz,t,beta1,beta2,n,tau,lambda,cov1,cov2))-as.numeric((Real_ExpDeriv1(Nlabel,phi,Di,EstAlpha,delta,lambda,i,zzz,t,EstBeta1,EstBeta2,n,tau,cov1,cov2))^2*Real_deriv2(Nlabel,EstU,Di,alpha1,delta,i,zzz,t,beta1,beta2,n,tau,lambda,cov1,cov2))
            }
          }
        }
      }
      I12[t]=sum(I11)
    }
    InfeNew12=sum(I12,na.rm=T)
    EndNew12=SusNew10+InfeNew12
    Estdelta=delta-EndNew10/EndNew12
    ##################################### Sigma_U #####################
    B=(lambda1*D+(1-lambda1)*I)
    firstS=zz/sigma1-sigma1*Real_Estsigma(EstU,D,lambda1,I)
    firstL=tr((D-I)%*%solve(B))/2+sigma1^2*Real_Estlambda(EstU,D,I)/2
    firstLS=sigma1*Real_Estlambda(EstU,D,I)
    secondS=-zz/sigma1^2-Real_Estsigma(EstU,D,lambda1,I)
    secondL=-tr((D-I)%*%solve(B)%*%(D-I)%*%solve(B))/2
    Inf1=matrix(c(secondS,firstLS,firstLS,secondL),2,2)
    Score=c(firstS,firstL)
    c(sigma1,lambda1)-solve(Inf1)%*%Score
    Aest=c(sigma1,lambda1)-solve(Inf1)%*%Score
    HatSigmmaU=Aest[1]
    EstGammau=Aest[2]
    ##############################################
    result=list(beta1=EstBeta1,beta2=EstBeta2,Uhat=EstU,alpha1=EstAlpha,
                delta=Estdelta,sigma1=HatSigmmaU,lambda1=EstGammau)
    result
  }

  RealData_Loglik=function(Nlabel,phi,Di,alpha1,delta,lambda1,sigma1,beta1,beta2,time,n,zz,tau,lambda,I,D,cov1,cov2){
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
              prob1=1-exp(-exp(alpha1+cov1[i,]%*%beta1+cov2[zzz,]
                               %*%beta2+phi[zzz])*dx)
              FL2[i]=log(1-prob1)
              FL2[i]=ifelse(is.infinite(FL2[i]),-20,FL2[i])
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
              prob1=1-exp(-exp(alpha1+cov1[i,]%*%beta1+cov2[zzz,]
                               %*%beta2+phi[zzz])*dx5)
              FL2[i]=log(prob1)
              FL2[i]=ifelse(is.infinite(FL2[i]),-20,FL2[i])
            }}}}
      FL1[t]=sum(FL2,na.rm=T)
    }
    Sigma1 = solve(sigma1^2*(lambda1*D+(1-lambda1)*I))
    LLL=dmvnorm(phi, rep(0,zz), solve(Sigma1*(as.numeric(lambda1)*D+(1-as.numeric(lambda1))*I)), log = FALSE)
    LLL1=ifelse(LLL>0,log(LLL),0)
    loglik1=(sum(FL1,na.rm=T))+LLL1
  }


  for (iter in 1:ITER) {
    mu = rep(0,zz)
    Sigma0=solve(sigma0^2*(lambda0*D+(1-lambda0)*I))
    phi1[iter,] = mvrnorm(1,mu,Sigma0, tol = 1e-6)
    phi = phi1[iter,]
    OUT2 = list()
    ss = 1
    LA = numeric()
    est0 = RealData_Estpar(Nlabel,phi,Di,alpha0,delta0,lambda0,sigma0,beta10,beta20,zz,time,n,tau,lambda,D,I,cov1,cov2,q1,q2)
    alpha1 = est0$alpha1
    alpha1
    delta = est0$delta
    delta
    lambda1 = est0$lambda1
    lambda1
    sigma1 = est0$sigma1
    sigma1
    beta1 = est0$beta1
    beta1
    beta2 = est0$beta2
    beta2
    Uhat = est0$Uhat
    LA[1] = RealData_Loglik(Nlabel,phi,Di,alpha1,delta,lambda1,sigma1,beta1,beta2,time,n,zz,tau,lambda,I,D,cov1,cov2)
    print(LA[1])
    OUT2[[1]] = c(alpha1, beta1, beta2, delta, sigma1, lambda1)
    repeat{
      ss = ss + 1
      if(ss > 16){
        alpha1 = NA
        delta = NA
        sigma1 = NA
        lambda1 = NA
        out1 = list(alpha1 = NA, delta = NA, sigma1 = NA, lambda1 = NA, beta1 = NA, beta2 = NA, ss = NA)
        break
      }
      if(sigma1 < 0 | lambda1 > 1 | lambda1 < 0){
        alpha1 = NA
        delta = NA
        sigma1 = NA
        lambda1 = NA
        out1 = list(alpha1 = NA, delta = NA, sigma1 = NA, lambda1 = NA, beta1 = NA, beta2 = NA, ss = ss)
        break
      }
      E0 = c(alpha1, beta1, beta2, delta, sigma1, lambda1)
      est = RealData_Estpar(Nlabel, phi, Di, alpha1, delta, lambda1, sigma1, beta1, beta2,zz,time,n,tau,lambda,D,I,cov1,cov2,q1,q2)
      alpha1 = est$alpha1
      beta1 = est$beta1
      beta2 = est$beta2
      delta = est$delta
      lambda1 = est$lambda1
      sigma1 = est$sigma1
      Uhat = est$Uhat
      if(sigma1 < 0 | lambda1 > 1 | lambda1 < 0){
        alpha1 = NA
        delta = NA
        sigma1 = NA
        lambda1 = NA
        out1 = list(alpha1 = NA, delta = NA, sigma1 = NA, lambda1 = NA, beta1 = NA, beta2 = NA, ss = ss, LA = LA)
        break
      }
      if(sigma1 >= 0 & lambda1 <= 1 & lambda1 >= 0){
        E1 = c(alpha1, beta1, beta2, delta, sigma1, lambda1)
        Dist = sqrt(sum((E0 - E1)^2))
        print(Dist)
        LA[ss] = RealData_Loglik(Nlabel,phi,Di,alpha1,delta,lambda1,sigma1,beta1,beta2,time,n,zz,tau,lambda,I,D,cov1,cov2)
        out1 = list(Uhat = Uhat, alpha1 = alpha1, beta1 = beta1, beta2 = beta2, delta = delta, sigma1 = sigma1, lambda1 = lambda1, ss = ss, LA = LA)
        OUT2[[ss]] = c(alpha1, beta1, beta2, delta, sigma1, lambda1)
        if(!abs((LA[ss] - LA[ss - 1])/LA[ss]) > 0.25){
          break
        }
      }
    }
    LoglikObs=function(Nlabel,phi,Di,alpha1,delta,beta1,beta2){
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
                prob1=1-exp(-exp(alpha1+cov1[i,]%*%beta1+cov2[zzz,]
                                 %*%beta2+phi[zzz])*dx)
                FL2[i]=log(1-prob1)
                FL2[i]=ifelse(is.infinite(FL2[i]),-20,FL2[i])
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
                prob1=1-exp(-exp(alpha1+cov1[i,]%*%beta1+cov2[zzz,]
                                 %*%beta2+phi[zzz])*dx5)
                FL2[i]=log(prob1)
                FL2[i]=ifelse(is.infinite(FL2[i]),-20,FL2[i])
              }}}}
        FL1[t]=sum(FL2,na.rm=T)
      }
      loglik1Obs=(sum(FL1,na.rm=T))
      return(loglik1Obs)
    }
    VecEM[iter,] = c(out1$alpha1, out1$beta1, out1$beta2, out1$delta, out1$sigma1, out1$lambda1)
    Finalphi[[iter]] = phi
    if(is.na(delta) && is.na(lambda1) && is.na(sigma1) && is.na(beta1) && is.na(beta2)){
      EstLoglikObs[iter] = NA
      RealLog[iter] = NA
      EstUhat[iter] = NA
      DiffValue[iter] = NA
    }
    if(!is.na(delta) && !is.na(lambda1) && !is.na(sigma1) ){
      EstLoglikObs[iter] = LoglikObs(Nlabel, phi, Di, out1$alpha1, out1$delta, out1$beta1, out1$beta2)
      RealLog[iter] = LoglikObs(Nlabel, Finalphi[[1]], Di, out1$alpha1, out1$delta, out1$beta1, out1$beta2)
      EstUhat[iter] = LoglikObs(Nlabel, Uhat, Di, out1$alpha1, out1$delta, out1$beta1, out1$beta2)
      DiffValue[iter] = EstUhat[iter] - RealLog[iter]
    }
  }

  fr=function(par){
    #par = c(alpha1,beta1,beta2,delta,sigma1,lambda1)
    alpha1=par[1]
    beta1[1,]=par[2]
    beta1[2,]=par[3]
    beta1[3,]=par[4]
    beta1[4,]=par[5]
    beta1[5,]=par[6]
    beta1[6,]=par[7]
    beta2[1,]=par[8]
    beta2[2,]=par[9]
    beta2[3,]=par[10]
    beta2[4,]=par[11]
    beta2[5,]=par[12]
    delta=par[13]
    sigma1=par[14]
    lambda1=par[15]

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
                  }
                }
              }
              dx1=replace(dx1,is.infinite(dx1),0)
              dx=sum(dx1)
              prob1=1-exp(-exp(alpha1+cov1[i,]%*%beta1+cov2[zzz,]
                               %*%beta2+Finalphi[[1]][zzz])*dx)
              FL2[i]=log(1-prob1)
              FL2[i]=ifelse(is.infinite(FL2[i]),-20,FL2[i])
            }
            if(tau[i]==(t+1)){
              dx4=rep(0,n)
              for(j in 1:n){
                if (Nlabel[j]==zzz){
                  if(tau[j]<=t & (tau[j]+lambda[j])>t & tau[j]!=0){
                    dx4[j]=Di[i,j]^(-delta)
                  }
                }
              }
              dx4=replace(dx4,is.infinite(dx4),0)
              dx5=sum(dx4)
              prob1=1-exp(-exp(alpha1+cov1[i,]%*%beta1+cov2[zzz,]
                               %*%beta2+Finalphi[[1]][zzz])*dx5)
              FL2[i]=log(prob1)
              FL2[i]=ifelse(is.infinite(FL2[i]),-20,FL2[i])
            }
          }
        }
      }
      FL1[t]=sum(FL2,na.rm=T)
    }
    Sigma1 = solve(sigma1^2*(lambda1*D+(1-lambda1)*I))
    LLL=dmvnorm(phi, rep(0,zz),
                solve(Sigma1^2*(as.numeric(lambda1)*D+(1-as.numeric(lambda1))*I)),
                log = TRUE)
    loglik1=((sum(FL1,na.rm=T))+LLL)
  }
  Est=VecEM[1,]
  x00=round(Est,2)
  hess <- hessian(fr,x00 )
  SE=sqrt(diag(solve(hess)))
  Est/SE
  ####################################Rate############################
  ###
  Rate1=array(0,c(n,time,zz))
  for(t in 1:time){
    for(zzz in 1:zz){
      for(i in 1:n){
        if (Nlabel[i]==zzz){
          if(tau[i]==(t+1)){
            dx1=rep(0,n)
            for(j in 1:n){
              if (Nlabel[j]==zzz){
                if(tau[j]<=t & (tau[j]+lambda[j])>t & tau[j]!=0){
                  dx1[j]=Di[i,j]^(-out1$delta)
                }
              }
            }
            dx2=replace(dx1,is.infinite(dx1),0)
            dx=sum(dx2,na.rm=T)
            Rate1[i,t,zzz]=1-exp(-(exp(Est[1]+cov1[i,]%*%Est[2:1+q1]+cov2[zzz,]
                                       %*%Est[q1+1:q1+q2+1]+Finalphi[[1]][zzz])*dx))
          }
        }
      }
    }
  }
  FinalRate1=c()


  }
