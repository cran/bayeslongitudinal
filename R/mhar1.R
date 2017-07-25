#' mhar1
#'
#' Run Bayesian estimation of a balanced longitudinal model with AR(1) structure
#' @param Data A vector with the observations of the response variable
#' @param Matriz The model design matrix
#' @param individuos A numerical value indicating the number of individuals in the study
#' @param tiempos A numerical value indicating the number of times observations were repeated
#' @param betai A vector with the initial values of the vector of regressors
#' @param rhoi A numerical value with the initial value of the correlation
#' @param beta1i A numerical value with the shape parameter of a beta apriori distribution of rho
#' @param beta2i A numerical value with the scaling parameter of a beta apriori distribution of rho
#' @param iteraciones A numerical value with the number of iterations that will be applied the algorithm MCMC
#' @param burn Number of iterations that are discarded from the chain
#' @return A dataframe with the mean, median and standard deviation of each parameter, A graph with the histograms and chains for the parameters that make up the variance matrix, as well as the selection criteria AIC, BIC and DIC
#' @references Gamerman, D. 1997. Sampling from the posterior distribution in generalized linear mixed models. Statistics and Computing, 7, 57-68
#' @references Cepeda, C and Gamerman, D. 2004. Bayesian modeling of joint regressions for the mean and covariance matrix. Biometrical journal, 46, 430-440.
#' @references  Cepeda, C and Nuñez, A. 2007. Bayesian joint modelling of the mean and covariance structures for normal longitudinal data. SORT. 31, 181-200.
#' @references  Nuñez A. and Zimmerman D. 2001. Modelación de datos longitudinales con estructuras de covarianza no estacionarias: Modelo de coeficientes aleatorios frente a modelos alternativos. Questio. 2001. 25.
#' @examples
#' attach(Dental)
#' Y=as.vector(distance)
#' X=as.matrix(cbind(1,age))
#' mhar1(Y,X,27,4,c(1,1),0.5,1,1,500,50)
#' @importFrom graphics hist par  plot
#' @importFrom stats dbeta median runif sd
#' @export
mhar1=function(Data,Matriz,individuos,tiempos,betai,rhoi,beta1i,beta2i,iteraciones,burn)
{
  Y=Data
  X=Matriz
  ind=individuos #need
  tie=tiempos #need
  Sigma=bloques(1,rhoi,tie,ind)#need
  nb=tie*ind
  contador=0
  iter=iteraciones
  ######Construccion de aprioris #################
  b0=betai
  B0=diag(100^4,nrow=length(b0),ncol=length(b0))
  g0=0.0001
  t0=0.0001
  s0=100
  rho_0=0.5
  beta1=beta1i  #need
  beta2=beta2i #need
  Mega=matrix(0,nrow=iter,ncol=length(b0))
  betas=matrix(0,nrow=iter,ncol=length(b0))
  medias=matrix(0,nrow=iter,ncol=nb)
  sigmas=rep(0,iter)
  prop=rep(0,iter)
  rho=rep(0,iter)
  dev=rep(0,iter)
  C=Sigma/s0
  #### Gibs media y varianza ####
  for(i in 2:iter){
    Mega[1,]=rep(0,length(b0))
    sigmas[1]=s0
    rho[1]=rho_0
    Sigma=C*sigmas[i-1]
    Bn=solve(t(X)%*%solve(Sigma)%*%X+solve(B0))
    Mega[i,]=Bn%*%(t(X)%*%solve(Sigma)%*%Y+solve(B0)%*%Mega[i-1,])
    betas[i,]=mvrnorm(1,Mega[i,],Bn)
    medias[i,]=as.vector(X%*%betas[i,])
    #####para varianza########
    Ra=t(Y-(X%*%betas[i,]))%*%(solve(C))%*%(Y-(X%*%betas[i,]))
    N=length(Y)
    shape=(N+g0)/2
    ratio=((g0*s0)+Ra)/2
    sigmas[i]=rigamma(1,shape,ratio)
    #shape=(N+g0)/2
    #ratio=(t0+Ra)/2
    #sigmas[i]=rigamma(1,shape,ratio)
    ####M-H####################
    #####Cadena############
    a=2*rho[i-1]
    if(rho[i-1]<=0.5)(prop[i]=runif(1,0,a))
    if(rho[i-1]>0.5)(prop[i]=runif(1,a-1,1))
    ############Bloques para C nuveo ###############################
    Cn=bloques(1,prop[i],tie,ind)
    ####################Bloques para C viejo ##################
    Cv=bloques(1,rho[i-1],tie,ind)
    ##############################################################
    #################Verosimilitud#######################
    Sigmav=Cv*sigmas[i]
    Sigman=Cn*sigmas[i]
    Lv=dmvnorm(as.vector(Y),medias[i,],Sigmav)
    Ln=dmvnorm(as.vector(Y),medias[i,],Sigman)
    pv=dbeta(rho[i-1],beta1,beta2)
    pn=dbeta(prop[i],beta1,beta2)
    acep=min(1,exp(log(Ln)+log(pn)-log(Lv)-log(pv)))
    if(runif(1)< acep )
    {rho[i]=prop[i]
    C=Cn
    contador=contador+1
    }else{
      rho[i]=rho[i-1]
      C=Cv
    }
    Sigma=C*sigmas[i]
    dev[i]=-2*log(dmvnorm(as.vector(Y),medias[i,],Sigma))
  }####for

  burnIn=burn

  ###############Para mostrar##############
  resulmedias=rep(0,(length(b0)+2))
  resulmedians=rep(0,(length(b0)+2))
  resulsd=rep(0,(length(b0)+2))
  for(k in 1:length(b0)){
    resulmedias[k]=mean(betas[(burnIn:iter),k])
    resulmedians[k]=median(betas[(burnIn:iter),k])
    resulsd[k]=sd(betas[(burnIn:iter),k])
  }
  resulmedias[length(b0)+1]=mean(sigmas[burnIn:iter])
  resulmedians[length(b0)+1]=median(sigmas[burnIn:iter])
  resulsd[length(b0)+1]=sd(sigmas[burnIn:iter])
  resulmedias[length(b0)+2]=mean(rho[burnIn:iter])
  resulmedians[length(b0)+2]=median(rho[burnIn:iter])
  resulsd[length(b0)+2]=sd(rho[burnIn:iter])
  Media1=as.vector(X%*%resulmedias[1:length(b0)])
  VARCOV=bloques(sqrt(mean(sigmas[burnIn:iter])),mean(rho[burnIn:iter]),tie,ind )
  LIK=dmvnorm(as.vector(Y),Media1,VARCOV)
  AIC=2*(length(b0)+2)-2*log(LIK)
  BIC=(log(nb)*(length(b0)+2))-2*log(LIK)
  #####para dic
  devm=-2*log(LIK)
  de=mean(dev[burnIn:iter])-devm
  DIC=devm+2*de
  tabla=data.frame(resulmedias,resulmedians,resulsd)
  letraa=0:(length(b0)-1)
  letrab=rep("b",length(b0))
  letraA=paste(letrab,letraa)
  letraA2=rep(0,length(letraA)+2)
  for(i in 1:length(letraA)){
    letraA2[i]=letraA[i]
  }
  letraA2[length(letraA)+1]="Var"
  letraA2[length(letraA)+2]="Cor"
  dimnames(tabla)=list(letraA2 , c("Mean","Median","S.D"))
  #return(tabla)
  par(mfrow=c(2,2))
  hist(sigmas[burnIn:iter],breaks=80,main=expression(sigma^2),xlab="",ylab="")
  plot(sigmas[burnIn:iter],type="l",xlab="",ylab="",main=expression(sigma^2))
  hist(rho[burnIn:iter],breaks=80,main=expression(rho),xlab="",ylab="")
  plot(rho[burnIn:iter],type="l",xlab="",ylab="",main=expression(rho))

  return(list(tabla,"AIC"=AIC, "BIC"=BIC,"DIC"=DIC))
}##function
