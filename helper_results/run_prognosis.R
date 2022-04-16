library(glmnet)
library(parallel)
library(survival)
library(data.table)
#library(mltools)
# library(CoxBoost)
library(randomForestSRC)
# library(CoxHD)
library(survC1)
### All the predictors
###---------------------------------------------------------------------------------------------------------------------------------------------------------------------

predictorGLM <- function(designTrain, designTest, responseTrain, alpha=0, ninternalfolds=10) {
    # alpha=1 --> l1 penalty
    # alpha=0 --> l2 penalty
    # alpha=1/2 --> elastic net
    set.seed(17)
    # Train
    cvfit = cv.glmnet(designTrain, responseTrain, family="cox", alpha=alpha, nfolds=ninternalfolds, grouped=TRUE)
    # Predict
    risk.predict = predict(cvfit, newx=designTest, s="lambda.min", type="response")
    risk.predict = as.vector(risk.predict[,1])

    return(risk.predict)
}    
predictorBoost<-function(designTrain, designTest, responseTrain, maxstepno=500 , K=10 , type="verweij" , penalty=100){
    
    set.seed(17)
    
    cv.res<-cv.CoxBoost(time=as.numeric(unlist(responseTrain[,1])),
                  status=as.numeric(unlist(responseTrain[,2])),
                  x=designTrain, maxstepno=maxstepno, K=K, type=type, penalty=penalty )
    
    cvfit <- CoxBoost(time=as.numeric(unlist(responseTrain[,1])),
                  status=as.numeric(unlist(responseTrain[,2])),
                  x=designTrain,stepno=cv.res$optimal.step , penalty=penalty)

    risk.predict <- predict(cvfit,designTest)
  
    return(as.vector(risk.predict))
}
predictorRF <- function(designTrain, designTest, responseTrain, ntree=100, importance="none",nodesize=10) {
    set.seed(17)
    # Train
    cvfit = rfsrc(Surv(time, status) ~ ., data=data.frame(designTrain,responseTrain), ntree=ntree, importance=importance,nodesize=nodesize)
    
    # Predict
    risk.predict = predict(cvfit, data.frame(designTest), importance=importance)$predicted
    
    return(risk.predict)
} 
predictorAIC <- function(designTrain, designTest, responseTrain) {
    set.seed(17)
    # Train
    c <- coxph(Surv(time, status) ~ ., data=data.frame(designTrain,responseTrain))
    scopeStep <- as.formula(paste("Surv(time,status) ~", paste(colnames(designTrain), collapse="+")))
    cvfit<-step(c, scope=scopeStep, k = 2, trace=0)
    # Predict
    risk.predict = predict(cvfit, data.frame(designTest))
    
    return(risk.predict)
}
predictorRFX <- function(designTrain, designTest, responseTrain, max.iter = 500,tol=0.01,groups = rep(1, ncol(designTrain))) {
    set.seed(17)
    # Train
    cvfit = CoxRFX(data.frame(designTrain), Surv(time=responseTrain[,1],event =responseTrain[,2]) , max.iter =max.iter,tol=tol,groups=groups)
    cvfit$Z <- NULL
    # Predict
    risk.predict<-predict(cvfit,data.frame(designTest))
    
    return(risk.predict)
}
predictorCPSS <- function(designTrain, designTest, responseTrain, bootstrap.samples = 50, nlambda = 100, mc.cores = 1) {
    set.seed(17)
    # Train
    cvfit = CoxCPSSInteractions(data.frame(designTrain), Surv(time=responseTrain[,1],event =responseTrain[,2]), bootstrap.samples = bootstrap.samples, nlambda = nlambda, mc.cores = mc.cores)
    # Predict
    risk.predict = predict(cvfit, data.frame(designTest))
    
    return(risk.predict)
}


###---------------------------------------------------------------------------------------------------------------------------------------------------------------------
### Cross Validation Function
###---------------------------------------------------------------------------------------------------------------------------------------------------------------------
runCV <- function(mypredictor, response, design, nfolds=nfolds, nrepeats=nrepeats, seed=seed, mc.cores=mc.cores, ...) {
    # function that run "mypredictor" on a CV setting
    #
    # output a list of size the number of CV experiments (eg 50) (= nfolds x nrepeats)
    
    # "ref" contains the responses of the fold test set

    #  random number generator seed
    set.seed(seed)
    # Make folds
    n = nrow(design)
    folds <- list()
    for (i in seq(nrepeats)) {
        folds <- c(folds,split(sample(seq(n)), rep(1:nfolds, length = n)))
    }
    nexp = length(folds) # the total number CV of experiments

    # Parallel CV
    print("start CV")
    rescv = mclapply(seq(nexp),
                   FUN=function(iexp) {
                       cat(".")
                       vTrain = design[-folds[[iexp]],,drop=F]
                       vTest = design[folds[[iexp]],,drop=F]
                       lTrain = response[-folds[[iexp]],]
                       lTest = response[folds[[iexp]],]
                       # Train and Predcit
                       #predict.test = ifelse(use_alpha==TRUE,mypredictor(designTrain=vTrain, designTest=vTest, response=lTrain,alpha=alpha, ...),mypredictor(designTrain=vTrain, designTest=vTest, response=lTrain, ...))
                       predict.test = mypredictor(designTrain=vTrain, designTest=vTest, response=lTrain, ...)
                       
                       # Evaluate CI on the test
                       ci.test = suppressWarnings(survConcordance(Surv(time,status) ~ predict.test, as.data.frame(lTest)))
                       model <- cbind(lTest,predict.test)

                       return(Est.Cval(model, tau=100,nofit=T)$Dhat)
                       #return(as.vector(ci.test$concordance))
                   },
                   mc.cores=mc.cores
                   )

    return(unlist(rescv))

}
###---------------------------------------------------------------------------------------------------------------------------------------------------------------------
### Launch prognosis
###---------------------------------------------------------------------------------------------------------------------------------------------------------------------
launch_prognosis <- function(x,y,predictors,str_predictors,nfolds=5,nrepeats=5,mc.cores=25,seed=17,...){
    set.seed(seed)
    i<-1
    str <- 1
    colnames(y) = c("time","status")
    table_w_c_i <- c()
    for(predictor in predictors){
        if (identical(predictorRF,predictor)){
            print("RandomForest")
            for (ntree in l_ntree){
                for (nodes in l_nodesize){
                
                    tmp <- runCV(mypredictor=predictor,
                                  response=y, design=x,
                                  nfolds=nfolds, nrepeats=nrepeats, seed=233,ntree=ntree,nodesize=nodes, mc.cores=mc.cores)
                    table_w_c_i <- cbind(table_w_c_i,tmp)
                    colnames(table_w_c_i) [i] <-paste(str_predictors[str],paste(ntree,nodes,sep="_"),sep="_")
                    i <- i+1
                    }
                }
                
        }
        
        if (identical(predictorGLM,predictor)){
            print("GLM")
            for (alpha in l_alpha){
                print("GLM")
                tmp <- runCV(mypredictor=predictor,
                                  response=y, design=x,
                                  nfolds=nfolds, nrepeats=nrepeats, seed=233,alpha=alpha, mc.cores=mc.cores)
                table_w_c_i <- cbind(table_w_c_i,tmp)
                colnames(table_w_c_i) [i] <-paste(str_predictors[str],alpha,sep="_") 
                i <- i+1
            }
        }
        
        if (identical(predictorBoost,predictor)){
            print("Boost")
             tmp <- runCV(mypredictor=predictor,
                              response=y, design=x,
                              nfolds=nfolds, nrepeats=nrepeats, seed=233,mc.cores=mc.cores,maxstepno=500 , K=10 , type="verweij" , penalty=100)
             table_w_c_i <- cbind(table_w_c_i,tmp)
             colnames(table_w_c_i) [i] <- str_predictors[str]
            i <- i+1
        }
        
        if (identical(predictorRFX,predictor)){
            print("RFX")
            tmp <- runCV(mypredictor=predictor,
                            response=y, design=x,
                            nfolds=nfolds, nrepeats=nrepeats, seed=233,mc.cores=mc.cores)
            table_w_c_i <- cbind(table_w_c_i,tmp)
            colnames(table_w_c_i) [i] <-str_predictors[str] 
            i <- i+1
        }
        str <- str +1
    }
    
    return (table_w_c_i)
    
    }
###---------------------------------------------------------------------------------------------------------------------------------------------------------------------