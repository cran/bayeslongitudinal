#' bloques3 compound symmetry
#'
#' Build a block diagonal matrix with compound symmetry structure
#' @param s Numerical value indicating global standard deviation of the matrix
#' @param r Numerical value indicating correlation of individuals
#' @param t Numerical value indicating number of times when observations are repeated
#' @param n Numerical value indicating number of individuals
#' @return A diagonal block matrix with compound symmetry structure
#' @references  Nuñez A. and Zimmerman D. 2001. Modelación de datos longitudinales con estructuras de covarianza no estacionarias: Modelo de coeficientes aleatorios frente a modelos alternativos. Questio. 2001. 25.
#' @examples
#' bloques3(2,.5,10,2)
#' @export
bloques3=function(s,r,t,n){
A=matrix(0,nrow=t,ncol=t)
for(i in 1:t){
for(j in 1:t){
if(i==j){
A[i,j]=s^2}
else{
A[i,j]=s^2*r
}
}##for2##
}##for1##
nb=t*n
m=matrix(0,nrow=nb,ncol=nb)
a=NULL
b=NULL
c=NULL
for(i in 1:n){
a[i]=i*t-(t-1)
b[i]=i*t
m[a[i]:b[i],a[i]:b[i]]=A
}
return(m)
}
#################################################################
#' bloques ar 1
#'
#' Build a block diagonal matrix with structure AR(1)
#' @param s Numerical value indicating global standard deviation of the matrix
#' @param r Numerical value indicating correlation of individuals
#' @param t Numerical value indicating number of times when observations are repeated
#' @param n Numerical value indicating number of individuals
#' @return A diagonal block matrix with structure AR(1)
#' @references  Nuñez A. and Zimmerman D. 2001. Modelación de datos longitudinales con estructuras de covarianza no estacionarias: Modelo de coeficientes aleatorios frente a modelos alternativos. Questio. 2001. 25.
#' @examples
#' bloques(2,.5,10,2)
#' @export

bloques=function(s,r,t,n){
A=matrix(0,nrow=t,ncol=t)
for(i in 1:t){
for(j in 1:t){
A[i,j]=s^2*r^(abs(i-j))
}##for2##
}##for1##
nb=t*n
m=matrix(0,nrow=nb,ncol=nb)
a=NULL
b=NULL
c=NULL
for(i in 1:n){
a[i]=i*t-(t-1)
b[i]=i*t
m[a[i]:b[i],a[i]:b[i]]=A
}
return(m)
}
##################################################################
#' bloques arma (1,1)
#'
#' Build a block diagonal matrix with structure ARMA(1,1)
#' @param s Numerical value indicating global standard deviation of the matrix
#' @param r Numerical value indicating the first parameter rho correlation of individuals
#' @param g Numerical value indicating the second parameter phi correlation of individuals
#' @param t Numerical value indicating number of times when observations are repeated
#' @param n Numerical value indicating number of individuals
#' @return A diagonal block matrix with structure ARMA(1,1)
#' @references  Nuñez A. and Zimmerman D. 2001. Modelación de datos longitudinales con estructuras de covarianza no estacionarias: Modelo de coeficientes aleatorios frente a modelos alternativos. Questio. 2001. 25.
#' @examples
#' bloques2(2,.5,.8,10,2)
#' @export
bloques2=function(s,r,g,t,n){
A=matrix(0,nrow=t,ncol=t)
for(i in 1:t){
for(j in 1:t){
if(i==j){
A[i,j]=s^2*r^(abs(i-j))}
else{
A[i,j]=s^2*r^(abs(i-j)-1)*g
}
}##for2##
}##for1##
nb=t*n
m=matrix(0,nrow=nb,ncol=nb)
a=NULL
b=NULL
c=NULL
for(i in 1:n){
a[i]=i*t-(t-1)
b[i]=i*t
m[a[i]:b[i],a[i]:b[i]]=A
}
return(m)
}

