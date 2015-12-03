\name{MethyImpute-package}
\alias{MethyImpute-package}
\alias{MethyImpute}
\docType{package}
\title{
\packageTitle{MethyImpute}
}
\description{
\packageDescription{MethyImpute}
}
\details{
The DESCRIPTION file:
\packageDESCRIPTION{MethyImpute}
\packageIndices{MethyImpute}
}
\author{
\packageAuthor{MethyImpute}

Maintainer: \packageMaintainer{MethyImpute}
}
\keyword{ package }

\examples{
#Do Not Run
## load the necessary packages
library("truncnorm")
library("mice")

## Set the initial conditions
N <- 500 ## Total sample size (subjects)
m <- 1000  ## Total CpG sites
missing.prop = 0.3 ## missing rate
impute.times = 10  ## The number of imputations
b4.value = 0.02    ## The effect size

seed = 1
data = simulate.data(seed,N,m,missing.prop, b4.value)

Y = data$reference.Y
pheno = data$pheno
pheno = pheno[,2:5]
reference.index = data$reference.num
missing.index = data$test.num
missing.cov.name = c("cell.1","cell.2","cell.3")
complete.cov.name = "x1"
m=30
maxit = 5
max.refernce.methy = 30
defaultMethod = c("norm", "logreg", "polyreg", "polr")
#imp = methy.mice(Y, pheno, missing.index, reference.index, missing.cov.name, complete.cov.name, max.refernce.methy = 30, m = 30, maxit = 5, defaultMethod = c("norm", "logreg", "polyreg", "polr"))
}
