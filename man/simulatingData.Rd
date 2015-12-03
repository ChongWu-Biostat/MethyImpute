\name{simulatingData}
\alias{simulatingData}
\title{
Generate Simulated Data Set}
\description{
Generate Simulated Data Set, which can be used to compare type 1 error rate and power for different methods}
\usage{
simulatingData(seed, N, m = 1000, missing.prop, b4.value)
}
\arguments{
  \item{seed}{random seed}
  \item{N}{
The number of subjects}
  \item{m}{
The number of CpG sites, the default is 1,000}
  \item{missing.prop}{
The missing rate}
  \item{b4.value}{
The effect sieze}
}

\value{
Return the simulated data set. 
}

\examples{

## Set the initial conditions
N <- 500 ## Total sample size (subjects)
m <- 1000  ## Total CpG sites
missing.prop = 0.3 ## missing rate
impute.times = 10  ## The number of imputations
b4.value = 0.02    ## The effect size

seed = 1
data = simulatingData(seed,N,m,missing.prop, b4.value)

}
