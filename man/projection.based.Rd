\name{projection.based}
\alias{projection.based}
\title{
Projection Based Method}
\description{
A extension to the Houseman 2012 method with adaptive reference sites.}
\usage{
projection.based(reference.CpG, reference.pheno, infer.CpG, missing.name, missing.n = 200)
}
\arguments{
  \item{reference.CpG}{Validation methylation (CpGs x subjects)}
  \item{reference.pheno}{Validation pheno}
  \item{infer.CpG}{Inference methylation (CpGs x subjects)}
  \item{missing.name}{The missing covariates name}
  \item{missing.n}{The number of CpG sites you want to use to impute the missing data}
}
\details{We tailor the method proposed by Houseman et al.to a more general situation. Briefly, the Houseman algorithm identifies 100-300 CpG sites that discriminated cellular composition in sorted normal human cell populations (consisting of B cells, granulocytes, monocytes, NK, and CD4+ and CD8+ T cells). The method fits a linear model at each of these CpG sites using a reference dataset to estimate the coefficient for each cellular component. It then uses a matrix projection approach to map these estimated coefficients to the relative proportions of each cellular component in the samples without cell type composition information. This method can be applied to other phenotypes besides cell type composition, if a reference dataset is available with both DNA methylation data and phenotype of interest. In the case of missing data, we use the set of samples with complete information as the reference dataset, and apply the projection approach to impute the missing values. 
}
\value{
Return the imputed covariates matrix
}
\author{
Chong Wu}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}

\keyword{projection.based}
