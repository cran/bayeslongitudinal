% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mhar1.R
\name{mhar1}
\alias{mhar1}
\title{mhar1}
\usage{
mhar1(Data, Matriz, individuos, tiempos, betai, rhoi, beta1i, beta2i,
  iteraciones, burn)
}
\arguments{
\item{Data}{A vector with the observations of the response variable}

\item{Matriz}{The model design matrix}

\item{individuos}{A numerical value indicating the number of individuals in the study}

\item{tiempos}{A numerical value indicating the number of times observations were repeated}

\item{betai}{A vector with the initial values of the vector of regressors}

\item{rhoi}{A numerical value with the initial value of the correlation}

\item{beta1i}{A numerical value with the shape parameter of a beta apriori distribution of rho}

\item{beta2i}{A numerical value with the scaling parameter of a beta apriori distribution of rho}

\item{iteraciones}{A numerical value with the number of iterations that will be applied the algorithm MCMC}

\item{burn}{Number of iterations that are discarded from the chain}
}
\value{
A dataframe with the mean, median and standard deviation of each parameter, A graph with the histograms and chains for the parameters that make up the variance matrix, as well as the selection criteria AIC, BIC and DIC
}
\description{
Run Bayesian estimation of a balanced longitudinal model with AR(1) structure
}
\examples{
attach(Dental)
Y=as.vector(distance)
X=as.matrix(cbind(1,age))
mhar1(Y,X,27,4,c(1,1),0.5,1,1,500,50)
}
\references{
Gamerman, D. 1997. Sampling from the posterior distribution in generalized linear mixed models. Statistics and Computing, 7, 57-68

Cepeda, C and Gamerman, D. 2004. Bayesian modeling of joint regressions for the mean and covariance matrix. Biometrical journal, 46, 430-440.

Cepeda, C and Nuñez, A. 2007. Bayesian joint modelling of the mean and covariance structures for normal longitudinal data. SORT. 31, 181-200.

Nuñez A. and Zimmerman D. 2001. Modelación de datos longitudinales con estructuras de covarianza no estacionarias: Modelo de coeficientes aleatorios frente a modelos alternativos. Questio. 2001. 25.
}
