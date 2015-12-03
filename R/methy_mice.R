methy.mice<-function(Y,pheno,missing.index, reference.index, missing.cov.name,complete.cov.name,max.refernce.methy = 30, m=30,maxit = 5,defaultMethod = c("norm", "logreg", "polyreg", "polr"))
{  #######################################################################
    ## set the simualtion parameter and get the simulation data set    ####
    ## I used the same way as Houseman 2014 reference free             ####
    #######################################################################
    test.num = missing.index
    reference.num = reference.index
    our.testData_Assay =Y[,test.num] ## out.testData_Assay, we need predict the cell type proportion for this part
    our.testData_Pheno =pheno[test.num,]
    
    ## our.validationData_Assay; using this part to impute the distribution of the cell type proportion
    our.validationData_Assay =  Y[,reference.num]
    our.validationData_Pheno = pheno[reference.num,]
    our.validationData_Assay = as.data.frame(our.validationData_Assay)
    tem.data = cbind(our.validationData_Pheno,t(our.validationData_Assay))

    pred.out = matrix(NA,length(missing.cov.name),dim(tem.data)[2])
    for(i in 1:length(missing.cov.name)) {
        tem.for = ""
        for(j in 1:dim(our.validationData_Assay)[1]) {
            if(j ==1){
                tem.for=  paste(tem.for,rownames(our.validationData_Assay)[j],sep="")
            } else {
                tem.for=  paste(tem.for,rownames(our.validationData_Assay)[j],sep="+")
            }
        }
        ### select the imputing variables for cell.1
        
        null.model = paste(missing.cov.name[i],"~1",sep = "")
        null = lm(as.formula(null.model),data = tem.data)
        full.model = paste(missing.cov.name[i],"~",complete.cov.name,"+",tem.for,sep = "")
        full = lm(as.formula(full.model),data = tem.data)
        eos.res = step(null,k= 2, scope=list(lower=null, upper=full), direction="forward",trace = 0)
        eos.index1 = names(eos.res$coefficients)
        eos.index1 = eos.index1[-1]

        if(length(eos.index1)>30) {
            eos.index1=eos.index1[1:30]
        }
        eos.index = !is.na(match(colnames(tem.data),eos.index1))
        
        pred.out[i,] = eos.index
    }

    ########### mice test
    ### Y1 is the who CpG sites, the first 250 are methylated, the reamining 750 are unmethylated
    test.data = cbind(pheno,t(Y))
    
    ## maxit = 0, get the default settings.
    ini <- mice(test.data, m = 25,maxit= 0,seed = 23109)
    
    pred <- ini$pred
    pred[missing.cov.name,] = 0
    pred[missing.cov.name,] = pred.out
    pred[missing.cov.name,complete.cov.name] = 1
    #pred[3,3] = 0
    ## m : number of multiple imputations
    ## use norm method to impute the missing data
    imp <- mice(test.data, m = m,maxit= maxit,predictorMatrix = pred,defaultMethod = defaultMethod)
    imp
}


methy.mice2<-function(Y,pheno,missing.index, reference.index, missing.cov.name,complete.cov.name, reference.methy.index, m=30,maxit = 5,defaultMethod = c("norm", "logreg", "polyreg", "polr"))
{  #######################################################################
    ## set the simualtion parameter and get the simulation data set    ####
    ## I used the same way as Houseman 2014 reference free             ####
    #######################################################################
    test.num = missing.index
    reference.num = reference.index
    our.testData_Assay =Y[,test.num] ## out.testData_Assay, we need predict the cell type proportion for this part
    our.testData_Pheno =pheno[test.num,]
    
    ## our.validationData_Assay; using this part to impute the distribution of the cell type proportion
    our.validationData_Assay =  Y[,reference.num]
    our.validationData_Pheno = pheno[reference.num,]
    our.validationData_Assay = as.data.frame(our.validationData_Assay)
    tem.data = cbind(our.validationData_Pheno,t(our.validationData_Assay))
    
    pred.out = matrix(NA,length(missing.cov.name),dim(tem.data)[2])
    for(i in 1:length(missing.cov.name)) {
        tem.for = ""
        for(j in 1:sum(!is.na(reference.methy.index[i,]))) {
            if(j ==1){
                tem.for=  paste(tem.for,rownames(Y)[reference.methy.index[i,j]],sep="")
            } else {
                tem.for=  paste(tem.for,rownames(Y)[reference.methy.index[i,j]],sep="+")
            }
        }
        ### select the imputing variables for cell.1
        
        null.model = paste(missing.cov.name[i],"~1",sep = "")
        null = lm(as.formula(null.model),data = tem.data)
        full.model = paste(missing.cov.name[i],"~",complete.cov.name,"+",tem.for,sep = "")
        full = lm(as.formula(full.model),data = tem.data)
        eos.res = step(null,k= 2, scope=list(lower=null, upper=full), direction="forward",trace = 0)
        
        eos.index1 = names(eos.res$coefficients)
        eos.index = !is.na(match(colnames(tem.data),eos.index1))
        pred.out[i,] = eos.index
    }
    
    ########### mice test
    ### Y1 is the who CpG sites, the first 250 are methylated, the reamining 750 are unmethylated
    test.data = cbind(pheno,t(Y))
    
    ## maxit = 0, get the default settings.
    ini <- mice(test.data, m = 25,maxit= 0,seed = 23109)
    
    pred <- ini$pred
    pred[missing.cov.name,] = 0
    pred[missing.cov.name,] = pred.out
    pred[missing.cov.name,complete.cov.name] = 1
    #pred[3,3] = 0
    ## m : number of multiple imputations
    ## use norm method to impute the missing data
    imp <- mice(test.data, m = m,maxit= maxit,predictorMatrix = pred,defaultMethod = defaultMethod)
    imp
}





