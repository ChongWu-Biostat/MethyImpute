# MethyImpute
Impute the missing covariates via Methylation Data

## Install the Package
We test it on R 3.2.1 in Linux server and 3.2.2 in Windows and Mac. *Note that for windows user, we need install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) first. For Mac user, we need install [gfortran](https://cran.r-project.org/bin/macosx/tools/).*

```
library(devtools)
install_github("ChongWu-Biostat/MethyImpute") # install the MethyImpute packages
```

## Simulated the Methylation and Covariates

Note that we simulate four cell type composite and one continuous covariate we are interested in. The methylation data come from three models. The first 50 CpG sites are correlated with the cell type composite and the covariate we are interested in. The 51-250 CpG sites are only correlated with the cell type composite, which is used to validate the type 1 error rate. The 251-400 CpG sites are only correlated with the covariate we are interested in. The 401-1000 CpG sites are drawn from the null distribution. 
```
library(MethyImpute)
## Set the initial conditions
N <- 500 ## Total sample size (subjects)
m <- 1000  ## Total CpG sites
missing.prop = 0.3 ## missing rate
b4.value = 0.02    ## The effect size
seed = 1
data = simulatingData(seed,N,m,missing.prop, b4.value)
```

## MI Based Method
Once we have the data, we need pre-process the data to the certain form and then we can apply multiple imputation (MI) based method to impute the missing covariate. 
```
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
## Note that this function is kind of time consuming. But it is feasiable to apply it to ARIC data or other real data set. 
imp = methy.mice(Y,pheno,missing.index, reference.index, missing.cov.name,complete.cov.name,max.refernce.methy = 30, m=30,maxit = 5,defaultMethod = c("norm", "logreg", "polyreg", "polr"))
```
Then we can apply the standard analysis. For example, we can use linear regression to test the relation between the methylation value and the covariate we are interested in as follows.
```
imputed.b4 = matrix(NA,5,5)
  for ( i in 1:5) {
    fit <- with(imp, lm(as.formula(paste(paste("methy",i,sep=""),"~ cell.1 + cell.2 + cell.3 +x1",sep =""))))
      fit.summary = summary(pool(fit))
      imputed.b4[i,1] = fit.summary[5,1] - 1.96 * fit.summary[5,2]
      imputed.b4[i,2] = fit.summary[5,1]
      imputed.b4[i,3] = fit.summary[5,1] + 1.96 * fit.summary[5,2]
      imputed.b4[i,4] = fit.summary[5,5] # p-vlaue
      imputed.b4[i,5] = fit.summary[5,2]
  }
  
  colnames(imputed.b4) = c("Lower","Estimated","Uppper","pvalue","SE")

```
## Projection Based Method
We tailor the method proposed by Houseman et al. to a more general situation. Briefly, the Houseman algorithm identifies 100-300 CpG sites that discriminated cellular composition in sorted normal human cell populations (consisting of B cells, granulocytes, monocytes, NK, and CD4+ and CD8+ T cells). The method fits a linear model at each of these CpG sites using a reference dataset to estimate the coefficient for each cellular component. It then uses a matrix projection approach to map these estimated coefficients to the relative proportions of each cellular component in the samples without cell type composition information. This method can be applied to other phenotypes besides cell type composition, if a reference dataset is available with both DNA methylation data and phenotype of interest. In the case of missing data, we use the set of samples with complete information as the reference dataset, and apply the projection approach to impute the missing values. 
```
Y = data$reference.Y
pheno = data$pheno

#################
reference.N = floor(N*(1-missing.prop))
reference.num = c(1:reference.N)

reference.Y = Y[,reference.num]
reference.pheno =pheno[reference.num,]

test.num = c((reference.N+1):N)
true.pheno = pheno[test.num,]
pheno[test.num,c(1:4)] =NA

our.testData_Assay =Y[,test.num]
our.testData_Pheno =pheno[test.num,]
our.validationData_Assay =  Y[,reference.num]
our.validationData_Pheno = pheno[reference.num,]


# reference.CpG contains the reference methylation data. We can run the single LMM to get the most significant CpG sites for each missing cell type proportion. Here we assume that we already have it.
reference.CpG =reference.Y  # Validation methylation (CpGs x subjects)

# The missing value name you want to impute
missing.name = c("cell.0","cell.1","cell.2","cell.3")


# number of CpG sites you want to use to impute the missing data
missing.n = 200

# The methylation data corresponding to the subjects that cell type proportion is missing
infer.CpG = our.testData_Assay # Inference methylation (CpGs x subjects)

##########################################################
## The main function to impute the missing data       ####
##########################################################
projection.based(reference.CpG,reference.pheno,infer.CpG,missing.name,missing.n = 200)
## NOTE: This function can directly apply to the other missing covariates, such as white blood cell count.
```


## Manual
If you like our pakcage *MethyImpute*, please give us a star. You can download the *MethyImpute* package manual [here](https://www.scribd.com/doc/292109233/MethyImpute-Manual). 
