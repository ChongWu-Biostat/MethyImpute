######################################################
# Multiple Imputation for cell type proportion    ####
# Version 1.00                                    ####
# Feb 2, 2015                                     ####
# Author: Chong Wu, Weihua Guan                   ####
######################################################


############################################################
### Simulation part Functions                           ####
############################################################
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


##############################################################################
##Define the simulation function                                         #####
## N: Total sample size                                                  #####
## m: Total CpG sites                                                    #####
## missing.prop: missing proportion                                      #####
##############################################################################
simulate.data <-function(seed,N,m=1000,missing.prop, b4.value)
{  #######################################################################
    ## set the simualtion parameter and get the simulation data set    ####
    ## I used the same way as Houseman 2014 reference free             ####
    #######################################################################
    
    #######################################################################
    simSeeds <- seed
    set.seed(simSeeds) ## Set the random seed for this simulation
    
    BETA.PROP <- 0.05
    BETA.NONZERO.MU <- 0.1
    BETA.NONZERO.SIG <- 0.025
    # Baseline proportion of cell types (e.g. T, B, mono, granulocyte)
    CELL.INTERCEPT <- c(0.2,0.1,0.1,0.6)  # Must sum to one
    
    # Cell mixture effect of slope
    CELL.SLOPE <- 0.5*c(0,-0.1,0.1,0) # must sum to zero
    
    CELL.DIRICHLETPRECISION <- 100
    
    MU.METHYLATION.BETA.A0 <- 0.25
    MU.METHYLATION.BETA.B0 <- 0.25
    N.METHYLATION.DMR <- 250
    
    EIG <- c(0.2,0.05)
    
    TOTALVARMICROARRAY <- 0.0625
    TOTALVAR <- 0.0325/4
    ERRDISPERSION <- 10
    
    N1 <- m   # Number of CpGs
    N2 <- N    # Number of subjects
    
    ZETA <- 1
    # Calculated quantities constant over all simulations
    NCELL <- length(CELL.SLOPE)  # Number of cell types
    
    ##
    b0 = rtruncnorm(N.METHYLATION.DMR, a=0, b=0.4, mean = 0.1, sd = 0.1)
    b1 = rtruncnorm(N.METHYLATION.DMR, a=0.5, b=0.95, mean = 0.7, sd = 0.1)
    b2 = rtruncnorm(N.METHYLATION.DMR, a=0.5, b=0.95, mean = 0.7, sd = 0.1)
    b3 = rtruncnorm(N.METHYLATION.DMR, a=0.5, b=0.95, mean = 0.7, sd = 0.1)
    
    MU.METHYLATION <- matrix( # matrix of cell-specific methylation mean values
    c(b0,b1,b2,b3),N.METHYLATION.DMR, NCELL)
    
    FIXEDEFF <- sampleDNAMethylCoeffs(N1,
    slope=list(mu=BETA.NONZERO.MU, sigma=BETA.NONZERO.SIG , pi=BETA.PROP),
    dmr=MU.METHYLATION%*%CELL.INTERCEPT)
    
    SIGMA <- TOTALVAR
    
    ALPHA.METH <- 1-1/(1+exp(FIXEDEFF[,1]))
    
    MU.METH <- outer(ALPHA.METH, rep(1,N2), "*")
    Fixed.Eff.x <- rep(0,N1)
    Fixed.Eff.x[1:(N.METHYLATION.DMR/5)] <-b4.value
    nonmethy.eff = (N1-N.METHYLATION.DMR)/5
    Fixed.Eff.x[(N.METHYLATION.DMR+1):(N.METHYLATION.DMR+nonmethy.eff)] <-b4.value
    BETA.METH <- ZETA*Fixed.Eff.x
    
    
    ## calculate the cell-type proportion
    x <- runif(N2, min=0.1, max=0.9)
    #x <- rtruncnorm(N2, a=0.1, b=0.9, mean = 0.5, sd = 0.3)
    
    # generate the cell type proportion
    dirichlet.wi = c(0.1,0.1,0.2,0.6)+as.matrix(c(0, 0,0.03,-0.03))%*% t(as.matrix(x))
    rou = 100
    omega = matrix(NA,N2,4)
    for (i in 1:N2)
    {
        omega[i,] = rdirichlet(1,(rou*dirichlet.wi)[,i])
    }
    ###############
    cell.0 <- omega[,1]
    cell.1 <- omega[,2]
    cell.2 <- omega[,3]
    cell.3 <- omega[,4]

    
    
    mu.cell <-  as.matrix(b0) %*% t(as.matrix(cell.0))+ as.matrix(b1) %*% t(as.matrix(cell.1))+as.matrix(b2) %*% t(as.matrix(cell.2))+as.matrix(b3) %*% t(as.matrix(cell.3))
    
    mu <- MU.METH
    # replace the methylation cell-type intercept with the effect of metylation cell type
    mu[1:N.METHYLATION.DMR,] <- mu.cell
    
    mu <- as.matrix(BETA.METH) %*% t(as.matrix(x)) + mu + (0.03*matrix(rnorm(N1*N2),N1,N2))
    
    #   mu <- ifelse(mu<0.000001, 0.000001, mu)
    #   mu <- ifelse(mu>0.999999, 0.999999, mu)
    
    muAvg <- apply(mu,1,mean)
    
    dispersion <- pmax(muAvg*(1-muAvg)/SIGMA/SIGMA,2)-1
    muA <- as.vector(dispersion*mu)
    muB <- as.vector(dispersion*(1-mu))
    
    # The Y is beta value, It's what we need.
    Y <- mu
    
    x1=x
    
    colnames(Y) =c(1:dim(Y)[2])
    rownames(Y) =c(1:dim(Y)[1])
    
    pheno = data.frame(cell.0, cell.1, cell.2, cell.3, x1)
    
    out.pheno = pheno[,2:4]
    
    ########################################################
    ### revised Multiple Imputation Method part        #####
    ########################################################
    ### Just the first 250 CpG sites are methylated
    rownames(Y) = paste("methy",1:dim(Y)[1],sep="")
    colnames(Y) = paste("sub",1:dim(Y)[2],sep="")
    rownames(pheno) = paste("sub",1:dim(pheno)[1],sep="")

    Y1 = Y   ## store the original data set
    Y = Y[1:250,]
    reference.N = floor(N*(1-missing.prop))  ## reference part means there is no missing data issue
    reference.num = c(1:reference.N)
    test.num = c((reference.N+1):N) ## test means the cell type proportion is missing
    true.pheno =   pheno[test.num,c(1:4)]
    whole.pheno = pheno
    

    pheno[test.num,c(1:4)] =NA
    
    test.data = cbind(pheno,t(Y))
    our.testData_Assay =Y[,test.num] ## out.testData_Assay, we need predict the cell type proportion for this part
    our.testData_Pheno =pheno[test.num,]
    
    ## our.validationData_Assay; using this part to impute the distribution of the cell type proportion
    our.validationData_Assay =  Y[,reference.num]
    our.validationData_Pheno = pheno[reference.num,]
    

    
    out = list(whole.Y = Y1,reference.Y = Y, pheno = pheno, whole.pheno = whole.pheno, reference.num = reference.num, test.num = test.num)
}



