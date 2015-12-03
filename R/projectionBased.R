###############################################
# Multiple Imputation for cell type proportion
#
#
# Version 1.10
# Oct 23, 2014
#
# Author: Chong Wu, Weihua Guan
#
################################################


##############################################################################
infer.interceptWBC <- function(Y, modelFix, coefWBC){
    
    
    M = dim(coefWBC)[1] #100
    n = dim(Y)[2]#184
    if(dim(Y)[1] != M) stop("Y does not match coefWBC in dimension\n")
    
    Z = cbind(1,coefWBC)
    Z = data.frame(Z)
    colnames(Z)[1] = "<Intercept>"
    Z$y = rep(0, M)
    
    #Y=coefEsts, X=Z
    xTest <- model.matrix(modelFix, Z)   #model.matrix creates a design (or model) matrix.
    sizeModel = dim(xTest)[2] #8
    
    # Initialize various containers
    nObserved = rep(NA, n)
    coefEsts = matrix(NA, n, sizeModel)
    coefVcovs =list()
    
    for(j in 1:n){ # For each CpG
        
        #Remove missing methylation values
        ii = !is.na(Y[,j])
        nObserved[j] = sum(ii)
        Z$y = Y[,j]
        
        if(j%%round(M/10)==0) cat(j,"\n") # Report progress
        # cat(j,"\n")
        # Try to fit a mixed model to adjust for plate
        fit = lm(modelFix, data=Z[ii,])
        coefEsts[j,] = fit$coef
        
    }
    # Name the rows so that they can be easily matched to the target data set
    rownames(coefEsts) = colnames(Y)
    colnames(coefEsts) = colnames(xTest)
    
    res = coefEsts
    sum =apply(res,1,sum)
    res2 =cbind(res,sum)
    return(res2)
}


validationWBC <- function(Y, pheno, modelFix, modelBatch=NULL, L.forFstat = NULL){
    N = dim(pheno)[1]
    pheno$y = rep(0, N)
    xTest <- model.matrix(modelFix, pheno)
    sizeModel = dim(xTest)[2]
    
    M = dim(Y)[1]
    
    if(is.null(L.forFstat)){
        L.forFstat = diag(sizeModel)
        colnames(L.forFstat) = colnames(xTest)
        rownames(L.forFstat) = colnames(xTest)
    }
    
    # Initialize various containers
    sigmaResid = sigmaIcept = nObserved = nClusters = Fstat = rep(NA, M)
    coefEsts = matrix(NA, M, sizeModel)
    coefVcovs =list()
    
    for(j in 1:M){ # For each CpG
        
        #Remove missing methylation values
        ii = !is.na(Y[j,])
        nObserved[j] = sum(ii)
        pheno$y = Y[j,]
        
        if(j%%round(M/10)==0) cat(j,"\n") # Report progress
        
        try({ # Try to fit a mixed model to adjust for plate
            if(!is.null(modelBatch)){
                fit = try(lme(modelFix, random=modelBatch, data=pheno[ii,]))
                OLS = inherits(fit,"try-error") # If LME can't be fit, just use OLS
            }
            else OLS = TRUE
            
            if(OLS){
                fit = lm(modelFix, data=pheno[ii,])
                fitCoef = fit$coef
                sigmaResid[j] = summary(fit)$sigma
                sigmaIcept[j] = 0
                nClusters[j] = 0
            }
            else{
                fitCoef = fit$coef$fixed
                sigmaResid[j] = fit$sigma
                sigmaIcept[j] = sqrt(getVarCov(fit)[1])
                nClusters[j] = length(fit$coef$random[[1]])
            }
            coefEsts[j,] = fitCoef
            coefVcovs[[j]] = vcov(fit)
            
            useCoef = L.forFstat %*% fitCoef
            useV = L.forFstat %*% coefVcovs[[j]] %*% t(L.forFstat)
            Fstat[j] = (t(useCoef) %*% solve(useV, useCoef))/sizeModel
        })
    }
    # Name the rows so that they can be easily matched to the target data set
    rownames(coefEsts) = rownames(Y)
    colnames(coefEsts) = names(fitCoef)
    
    degFree = nObserved - nClusters - sizeModel + 1
    
    # Get P values corresponding to F statistics
    Pval = 1-pf(Fstat, sizeModel, degFree)
    
    out = list(coefEsts=coefEsts, coefVcovs=coefVcovs, modelFix=modelFix, modelBatch=modelBatch,
    sigmaIcept=sigmaIcept, sigmaResid=sigmaResid, L.forFstat=L.forFstat, Pval=Pval, orderFstat=order(-Fstat),
    Fstat=Fstat, nClusters=nClusters, nObserved=nObserved, degFree=degFree)
    
    out
}

#########################################################

##########################################################
projection.based <-function(reference.CpG,reference.pheno,infer.CpG,missing.name,missing.n = 200)
{
    
    if(missing.n > dim(reference.CpG)[1]) {
        stop("The number of CpG sites for imputation is larger than the maximum number of CpG sites you provide.")
    }
    missing.name.fianl = ""
    
    for(i in 1:length(missing.name)) {
        if(i ==1){
            missing.name.fianl=  paste(missing.name.fianl,missing.name[1],sep="")
        } else {
            missing.name.fianl=  paste(missing.name.fianl,missing.name[i],sep="+")
        }
    }
    imputing.formula = as.formula(paste("y~",missing.name.fianl,"-1",sep =""))


    validEst = validationWBC(
    reference.CpG,             # Validation methylation (CpGs x subjects)
    reference.pheno,             # Validation phenotype frame (subjects x covariates)
    imputing.formula, # Validation model (fixed effects) [*Footnote 2*]
    #~1|BeadChip                       # Validation batch adjustment (random effects) [*Footnote 3*]
    )
    validEst
    
    # Choose top pseudo-DMRs
    DMRselection = rownames(validEst$coefEsts[validEst$orderFstat[1:missing.n],])
   
    OmegaEst = infer.interceptWBC(infer.CpG[DMRselection,],
    imputing.formula,
    validEst$coefEsts[DMRselection,])
    
    return(OmegaEst)
}
############################################


###########################################################################
## The following function is used for generating the simulated data     ###
###########################################################################
rdirichlet<-function(n,alpha)
## generate n random deviates from the Dirichlet function with shape
## parameters alpha
{
    l<-length(alpha);
    x<-matrix(rgamma(l*n,alpha),ncol=l,byrow=TRUE);
    sm<-x%*%rep(1,l);
    x/as.vector(sm);
}


sampleDNAMethylCoeffs <-function(
n,
intercept=list(
mu=c(-2.8, -0.75, 1.5),
sigma=c(0.7, 1.4, 0.7),
pi=c(0.3, 0.2, 0.5)
),
slope=list(
mu=0.4, sigma=0.1, pi=0.2
),
dmr=NULL
){
    comp <- sample(1:3, n, replace=TRUE, prob=intercept$pi)
    coef1 <- rnorm(n)*intercept$sigma[comp] + intercept$mu[comp]
    
    if(!is.null(dmr)) coef1[1:length(dmr)] <- dmr
    piEff <- slope$pi/intercept$pi[2]
    
    comp2 <- which(comp==2)
    ncomp2 <- length(comp2)
    coef2 <- rep(0, n)
    coef2[comp2] <- (-sign(coef1[comp2]-intercept$mu[2])*
    rbinom(ncomp2, 1, piEff)*rnorm(ncomp2, slope$mu, slope$sigma))
    cbind(coef1, coef2)
}

getRemainingMethylVariance <- function(d, sigma=0.3){
    d2sum <- sum(d*d)
    s2 = sigma*sigma
    s2-d2sum
}


