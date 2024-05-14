#' This function calculates the value of the log-likelihood function
#'
#' @param Nlabel Label for each sample from the area
#' @param phi Spatial random effects
#' @param Di Euclidean distance between susceptible individual and infectious individual
#' @param alpha1 intercept
#' @param delta Spatial parameter
#' @param lambda a vector containing the length of infectious period
#' @param sigma1 precision of spatial random effects
#' @param beta1 the parameter corresponding to the covariate associated with susceptible individual
#' @param beta2 the parameter corresponding to the covariate associated with area
#' @param time time
#' @param n Total number of individuals
#' @param zz Number of areas
#' @param tau the set of infectious individuals at time t in the zth area
#' @param lambda1 Spatial dependence
#' @param I Identity matrix
#' @param D matrix reflecting neighborhood structure
#' @param cov1 Individual level covariates
#' @param cov2 Area level covariates
#'
#' @return a numeric value for the log-likelihood
#' @export
#'
#'
#' @examples Sim_Loglik(rep(1:4,each=5), runif(4,min = 0, max = 1),
#' matrix(runif(400,min=4,max=20),nrow=20,byrow=TRUE),0.4, 2,rep(3,20),0.5,1,1,
#' 10,20,4,sample(c(0,1),replace = TRUE, size = 20),0.6,diag(4),
#' matrix(c(0,-1,0,-1,-1,0,-1,-1,0,-1,0,-1,-1,-1,-1,0),nrow=4,byrow=TRUE),
#' runif(20, 0, 1), runif(4, 0, 1))
#'
Sim_Loglik=function(Nlabel,phi,Di,alpha1,delta,lambda,sigma1,beta1,beta2,time,n,
                    zz,tau,lambda1,I,D,cov1,cov2){
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
  LLL=dmvnorm(phi, rep(0,zz), solve(Sigma1*(as.numeric(lambda1)*D+(1-as.numeric(lambda1))*I)), log = FALSE)
  LLL1=ifelse(LLL>0,log(LLL),0)
  loglik1=(sum(FL1,na.rm=T))+LLL1
}
