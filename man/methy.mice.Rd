\name{methy.mice}
\alias{methy.mice}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Multiple Imputation Based Method
}
\description{
		Multiple Imputation Based Method with two step variable selection 
}
\usage{
methy.mice(Y, pheno, missing.index, 
reference.index, missing.cov.name, complete.cov.name,
 max.refernce.methy = 30, m = 30, maxit = 5, 
 defaultMethod = c("norm", "logreg", "polyreg", "polr"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{Methylation Data Set}
  \item{pheno}{A data frame or a matrix containing the incomplete data.  Missing
values are coded as \code{NA}.}
  \item{missing.index}{Indicators of the subjects who have the missing covariates}
  \item{reference.index}{Indicators of the subjects who have the complete covariates}
  \item{missing.cov.name}{The missing covariate name}
  \item{complete.cov.name}{The complete covariate name}
  \item{max.refernce.methy}{The maximum CpG sites will be used to impute the missing covariates}
  \item{m}{Number of multiple imputations. The default is \code{m=5}.}
  \item{maxit}{A scalar giving the number of iterations. The default is 5.}
  \item{defaultMethod}{A vector of three strings containing the default
imputation methods for numerical columns, factor columns with 2 levels, and
columns with (unordered or ordered) factors with more than two levels,
respectively. If nothing is specified, the following defaults will be used:
\code{pmm}, predictive mean matching (numeric data) \code{logreg}, logistic
regression imputation (binary data, factor with 2 levels) \code{polyreg},
polytomous regression imputation for unordered categorical data (factor >= 2
levels) \code{polr}, proportional odds model for (ordered, >= 2 levels)}
}
\details{Our proposed method combines multiple imputation and variable selection in high-dimensional data. To deal with high-dimensional methylation data, we use a two-step variable selection approach including a screen and a selection stage. Then standard multiple imputation approach is applied to impute the missing covariate values and account for the uncertainty of imputation. 
	
Generates multiple imputations for incomplete multivariate data by Gibbs
sampling. Missing data can occur anywhere in the data. The algorithm imputes
an incomplete column (the target column) by generating 'plausible' synthetic
values given other columns in the data. Each incomplete column must act as a
target column, and has its own specific set of predictors. The default set of
predictors for a given target consists of all other columns in the data. For
predictors that are incomplete themselves, the most recently generated
imputations are used to complete the predictors prior to imputation of the
target column.

A separate univariate imputation model can be specified for each column. The
default imputation method depends on the measurement level of the target
column. In addition to these, several other methods are provided. You can
also write their own imputation functions, and call these from within the
algorithm.

The data may contain categorical variables that are used in a regressions on
other variables. The algorithm creates dummy variables for the categories of
these variables, and imputes these from the corresponding categorical
variable. The extended model containing the dummy variables is called the
padded model. Its structure is stored in the list component \code{pad}.

Built-in elementary imputation methods are:

\describe{
\item{pmm}{Predictive mean matching (any)}
\item{norm}{Bayesian linear regression (numeric)}
\item{norm.nob}{Linear regression ignoring model error (numeric)}
\item{norm.boot}{Linear regression using bootstrap (numeric)}
\item{norm.predict}{Linear regression, predicted values (numeric)}
\item{mean}{Unconditional mean imputation (numeric)}
\item{2l.norm}{Two-level normal imputation (numeric)}
\item{2l.pan}{Two-level normal imputation using pan (numeric)}
\item{2lonly.mean}{Imputation at level-2 of the class mean (numeric)}
\item{2lonly.norm}{Imputation at level-2 by Bayesian linear regression (numeric)}
\item{2lonly.pmm}{Imputation at level-2 by Predictive mean matching (any)}
\item{quadratic}{Imputation of quadratic terms (numeric)}
\item{logreg}{Logistic regression (factor, 2 levels)}
\item{logreg.boot}{Logistic regression with bootstrap}
\item{polyreg}{Polytomous logistic regression (factor, >= 2 levels)}
\item{polr}{Proportional odds model (ordered, >=2 levels)}
\item{lda}{Linear discriminant analysis (factor, >= 2 categories)}
\item{cart}{Classification and regression trees (any)}
\item{rf}{Random forest imputations (any)}
\item{ri}{Random indicator method for nonignorable data (numeric)}
\item{sample}{Random sample from the observed values (any)}
\item{fastpmm}{Experimental: Fast predictive mean matching using C++ (any)}
}
}
\author{
Chong Wu}
\examples{
### Do Not Run
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methy.mice}
