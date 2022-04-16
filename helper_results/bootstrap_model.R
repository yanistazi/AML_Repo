library(glmnet)
library(survC1)
library(survival)
options(warn=-1)




################################## A. NCRI Bootstrap comparison of feature combinations ##################################

# Prepare dataset

df_training <- read.table("aml_prognosis_updated.tsv")
normalized<-function(y) {

  x<-y[!is.na(y)]

  x<-(x - min(x)) / (max(x) - min(x))

  y[!is.na(y)]<-x

  return(y)
}

df_training$gender <- factor(df_training$gender,labels=c(0,1))
df_training$ahd <- factor(df_training$ahd,labels=c(0,1))
df_training$perf_status <- factor(df_training$perf_status,labels=c(0,1,2,3,4))
df_training$secondary <- factor(df_training$secondary,labels=c(1,2,3))
df_training[,c("age","bm_blasts","wbc","hb","plt")] <- apply(df_training[,c("age","bm_blasts","wbc","hb","plt")],2,normalized) 

eln <- c(2,3,4)
age <- c(167)

all_gen <- c(5:88)
vect <- apply(X=df_training[,all_gen],2,FUN=function(x) 100*length(which(x==1))/dim(df_training)[1])
gen <- match(names(vect[vect>=2]),names(df_training))
              
all_cyto <- c(89:158)
vect <- apply(X=df_training[,all_cyto],2,FUN=function(x) 100*length(which(x==1))/dim(df_training)[1])
cyto <- match(names(vect[vect>=2]),names(df_training))
              
clin <- c(159:165)
demo <- c(166:167)

comp <-c(190:205)                        

clin_demo <-c(clin,demo)
              
all_gen_cyto <- c(all_gen,all_cyto)
              
all_gen_cyto_clin_demo <- c(all_gen_cyto,clin_demo)

comp_ITD <- c(comp,29)
              
comp_ITD_clin_demo <- c(comp_ITD,clin,demo)
              
comp_clin_demo <- c(comp,clin,demo)
              
df_training$time <- df_training$os
df_training$status <- df_training$os_status

              
########### Split the data into test and validation (75-25)
# For reproducibility, we use the same seed
IDTRAIN <- readRDS("list_NCRI_Cohort_Training.rds") 
IDTEST <- readRDS("list_NCRI_Cohort_Testing.rds") 

dataTRAIN <- df_training[IDTRAIN,]
dataTEST <- df_training[IDTEST,]

              
              
#convert factor to numeric for the matrix multiplication
dataTEST$ahd <- as.numeric(as.character(dataTEST$ahd))
dataTEST$perf_status <- as.numeric(as.character(dataTEST$perf_status))
dataTEST$secondary <- as.numeric(as.character(dataTEST$secondary))
dataTEST$gender <- as.numeric(as.character(dataTEST$gender))


yTRAIN <- dataTRAIN[,c("time","status")]


prognosis_features<- list(eln=eln,all_gen=all_gen,comp=comp,comp_ITD=comp_ITD,
                          all_gen_cyto=all_gen_cyto,clin_demo=clin_demo,all_gen_cyto_clin_demo=all_gen_cyto_clin_demo,comp_ITD_clin_demo=comp_ITD_clin_demo)





coeffs <- list()
dataTESTmodels <- list() 
bscmodels <- list()

bootstrap_iter <- 100

for (i in 1:length(prognosis_features)){
    coeffs[names(prognosis_features)[i]] <- coef(cv.glmnet(data.matrix(dataTRAIN[,prognosis_features[[i]]]), data.matrix(yTRAIN), family="cox", alpha=0, nfolds=5, grouped=TRUE))
    RS <- data.matrix(dataTEST[,prognosis_features[[i]]]) %*% coeffs[[names(prognosis_features)[i]]]
    dataTESTmodels[names(prognosis_features)[i]] <- cbind(dataTEST$time, dataTEST$status, RS)
    bscmod <- NULL
    for(bs in 1:bootstrap_iter){
        set.seed(bs)  # so that all feature combination use same bootstrap data
        bssample <- sample(1:nrow(dataTEST), replace=TRUE)

        bscmod <- c(bscmod,Est.Cval(dataTESTmodels[[names(prognosis_features)[i]]][bssample,], tau=100,nofit=TRUE)$Dhat)
    }
    bscmodels[[names(prognosis_features)[i]]] <- bscmod
}

write.table(bscmodels,"bootstrap_results_NCRI_Cohort.tsv")


################################## B. NEJM Cohort Bootstrap comparison of feature combinations ##################################

library(glmnet)
library(survC1)
library(survival)
options(warn=-1)

df_validation <- read.table("full_data_validation.tsv")

normalized<-function(y) {

  x<-y[!is.na(y)]

  x<-(x - min(x)) / (max(x) - min(x))

  y[!is.na(y)]<-x

  return(y)
}


df_validation$Gender <- factor(df_validation$Gender,labels=c(0,1))
df_validation$perf_status <- factor(df_validation$perf_status,labels=c(0,1,2,3,4))
df_validation$ahd <- factor(df_validation$ahd,labels=c(0,1))
df_validation$secondary <- factor(df_validation$secondary,labels=c(1,2,3))
df_validation[,c("Age","BM_Blasts","WBC","HB","PLT")] <- apply(df_validation[,c("Age","BM_Blasts","WBC","HB","PLT")],2,normalized) 

eln <- c(111,112,113)


all_gen <- c(1:57)
vect <- apply(X=df_validation[,all_gen],2,FUN=function(x) 100*length(which(x==1))/dim(df_validation)[1])
gen <- match(names(vect[vect>=2]),names(df_validation))

all_cyto <- c(58:81)
vect <- apply(X=df_validation[,all_cyto],2,FUN=function(x) 100*length(which(x==1))/dim(df_validation)[1])
cyto <- match(names(vect[vect>=2]),names(df_validation))

clin <- c(84:90)
demo <- c(82:83)

comp <- c(115:130)
comp_ITD <- c(comp,23)

all_gen_cyto <- c(all_gen,all_cyto)
gen_cyto <- c(gen,cyto)
clin_demo <- c(clin,demo)
all_gen_cyto_clin_demo <- c(all_gen_cyto,clin_demo)
gen_cyto_clin_demo <- c(gen_cyto,clin_demo)
comp_ITD_clin_demo <- c(comp_ITD,clin_demo)


df_validation$time <- df_validation$OS
df_validation$status <- df_validation$OS_Status

########### Split the data into test and validation (75-25)
# For reproducibility, we use the same seed
              
IDTRAIN <- readRDS("list_NEJM_Cohort_Training.rds") 
IDTEST <- readRDS("list_NEJM_Cohort_Testing.rds") 

dataTRAIN <- df_validation[IDTRAIN,]
dataTEST <- df_validation[IDTEST,]

#convert factor to numeric for the matrix multiplication
dataTEST$ahd <- as.numeric(as.character(dataTEST$ahd))
dataTEST$perf_status <- as.numeric(as.character(dataTEST$perf_status))
dataTEST$secondary <- as.numeric(as.character(dataTEST$secondary))
dataTEST$Gender <- as.numeric(as.character(dataTEST$Gender))


yTRAIN <- dataTRAIN[,c("time","status")]




prognosis_features<- list(eln=eln,all_gen=all_gen,comp=comp,comp_ITD=comp_ITD,
                          all_gen_cyto=all_gen_cyto,clin_demo=clin_demo,all_gen_cyto_clin_demo=all_gen_cyto_clin_demo,comp_ITD_clin_demo=comp_ITD_clin_demo)


coeffs <- list()
dataTESTmodels <- list() 
bscmodels <- list()

bootstrap_iter <- 100

for (i in 1:length(prognosis_features)){
    coeffs[names(prognosis_features)[i]] <- coef(cv.glmnet(data.matrix(dataTRAIN[,prognosis_features[[i]]]), data.matrix(yTRAIN), family="cox", alpha=0, nfolds=5, grouped=TRUE))
    RS <- as.matrix(dataTEST[,prognosis_features[[i]]]) %*% coeffs[[names(prognosis_features)[i]]]
    dataTESTmodels[names(prognosis_features)[i]] <- cbind(dataTEST$time, dataTEST$status, RS)
    bscmod <- NULL
    for(bs in 1:bootstrap_iter){
        set.seed(bs)  # so that all feature combination use same bootstrap data
        bssample <- sample(1:nrow(dataTEST), replace=TRUE)
        bscmod <- c(bscmod,Est.Cval(dataTESTmodels[[names(prognosis_features)[i]]][bssample,], tau=100,nofit = TRUE)$Dhat)
    }
    bscmodels[[names(prognosis_features)[i]]] <- bscmod
}

write.table(bscmodels,"bootstrap_results_NEJM_Cohort.tsv")


################################## C. NCRI Comparison ELN vs Proposal ##################################

molecular_proposal_classification <- function(df){
    
    df$molecular_classification <- "none"

    df[(df$principal_component_NPM1==1 |df$principal_component_inv_16==1 |
                df$principal_component_t_8_21==1 | df$principal_component_t_15_17==1  |
                df$principal_component_CEBPA_bi==1 | df$principal_component_no_events==1 ),"molecular_classification"] <- "NEW_favorable"

    df[(df$principal_component_sAML1==1 | df$principal_component_t_6_9==1 | df$principal_component_DNMT3A_IDH1_2==1 |
               df$principal_component_Trisomies==1 | df$principal_component_t_11==1 | df$principal_component_WT1==1 |
                df$principal_component_DNMT3A_IDH1_2==1 | df$principal_component_mNOS==1),"molecular_classification"] <- "NEW_intermediate"

    df[( df$principal_component_sAML2==1 | df$principal_component_TP53_complex==1   |
                df$principal_component_inv_3==1) ,"molecular_classification"] <- "NEW_adverse"
    
    df[df$molecular_classification=="NEW_intermediate" & df$ITD==1,"molecular_classification"] <- "NEW_adverse"

    df[df$molecular_classification=="NEW_favorable" & df$principal_component_NPM1==1 & df$ITD==1,"molecular_classification"] <- "NEW_intermediate"

    
    return(df)

} 

library(glmnet)
library(survC1)
library(survival)
options(warn=-1)


#### Prepare proposal

df_training <- read.table("aml_prognosis_updated.tsv")
df_training <- molecular_proposal_classification(df_training)


df_training$NEW_favorable <-  ifelse(df_training$molecular_classification=="NEW_favorable",1,0)
df_training$NEW_intermediate <-  ifelse(df_training$molecular_classification=="NEW_intermediate",1,0)
df_training$NEW_adverse <-  ifelse(df_training$molecular_classification=="NEW_adverse",1,0)
table(df_training$molecular_classification)


df_training$time <- df_training$os
df_training$status <- df_training$os_status

########### Split the data into test and validation (75-25)
# For reproducibility, we use the same seed

IDTRAIN <- readRDS("list_NCRI_Cohort_Training.rds") 
IDTEST <- readRDS("list_NCRI_Cohort_Testing.rds") 

dataTRAIN <- df_training[IDTRAIN,]
dataTEST <- df_training[IDTEST,]




yTRAIN <- dataTRAIN[,c("time","status")]

molecular_classification <- c("NEW_favorable","NEW_intermediate","NEW_adverse")
eln_2017 <- c("eln_2017_favorable","eln_2017_intermediate","eln_2017_adverse")


prognosis_features<- list(molecular_classification=molecular_classification,eln_2017=eln_2017)


coeffs <- list()
dataTESTmodels <- list() 
bscmodels <- list()

bootstrap_iter <- 100

for (i in 1:length(prognosis_features)){
    coeffs[names(prognosis_features)[i]] <- coef(cv.glmnet(as.matrix(dataTRAIN[,prognosis_features[[i]]]), as.matrix(yTRAIN), family="cox", alpha=0, nfolds=5, grouped=TRUE))
    RS <- as.matrix(dataTEST[,prognosis_features[[i]]]) %*% coeffs[[names(prognosis_features)[i]]]
    dataTESTmodels[names(prognosis_features)[i]] <- cbind(dataTEST$time, dataTEST$status, RS)
    bscmod <- NULL
    bscrsq <- NULL
    for(bs in 1:bootstrap_iter){
        set.seed(bs)  # so that all feature combination use same bootstrap data
        bssample <- sample(1:nrow(dataTEST), replace=TRUE)

        bscmod <- c(bscmod,Est.Cval(dataTESTmodels[[names(prognosis_features)[i]]][bssample,], tau=100)$Dhat)
        

    }
    bscmodels[[names(prognosis_features)[i]]] <- bscmod

}

write.table(bscmodels,"bootstrap_results_NCRI_cohort_eln_proposal.tsv")


################################## D. NEJM Comparison ELN vs Proposal ##################################

molecular_proposal_classification <- function(df){
    
    df$molecular_classification <- "none"

    df[(df$principal_component_NPM1==1 |df$principal_component_inv_16==1 |
                df$principal_component_t_8_21==1 | df$principal_component_t_15_17==1  |
                df$principal_component_CEBPA_bi==1 | df$principal_component_no_events==1 ),"molecular_classification"] <- "NEW_favorable"

    df[(df$principal_component_sAML1==1 | df$principal_component_t_6_9==1 | df$principal_component_DNMT3A_IDH1_2==1 |
               df$principal_component_Trisomies==1 | df$principal_component_t_11==1 | df$principal_component_WT1==1 |
                df$principal_component_DNMT3A_IDH1_2==1 | df$principal_component_mNOS==1),"molecular_classification"] <- "NEW_intermediate"

    df[( df$principal_component_sAML2==1 | df$principal_component_TP53_complex==1   |
                df$principal_component_inv_3==1) ,"molecular_classification"] <- "NEW_adverse"
    
    df[df$molecular_classification=="NEW_intermediate" & df$ITD==1,"molecular_classification"] <- "NEW_adverse"

    df[df$molecular_classification=="NEW_favorable" & df$principal_component_NPM1==1 & df$ITD==1,"molecular_classification"] <- "NEW_intermediate"

    
    return(df)

}  

library(glmnet)
library(survC1)
library(survival)
options(warn=-1)


#### Prepare proposal

df_validation <- read.table("full_data_validation.tsv")
df_validation <- molecular_proposal_classification(df_validation)


df_validation$NEW_favorable <-  ifelse(df_validation$molecular_classification=="NEW_favorable",1,0)
df_validation$NEW_intermediate <-  ifelse(df_validation$molecular_classification=="NEW_intermediate",1,0)
df_validation$NEW_adverse <-  ifelse(df_validation$molecular_classification=="NEW_adverse",1,0)
table(df_validation$molecular_classification)


df_validation$time <- df_validation$OS
df_validation$status <- df_validation$OS_Status

########### Split the data into test and validation (75-25)
# For reproducibility, we use the same seed

IDTRAIN <- readRDS("list_NEJM_Cohort_Training.rds") 
IDTEST <- readRDS("list_NEJM_Cohort_Testing.rds") 

dataTRAIN <- df_validation[IDTRAIN,]
dataTEST <- df_validation[IDTEST,]




yTRAIN <- dataTRAIN[,c("time","status")]

molecular_classification <- c("NEW_favorable","NEW_intermediate","NEW_adverse")
eln_2017 <- c("eln_2017_favorable","eln_2017_intermediate","eln_2017_adverse")


prognosis_features<- list(molecular_classification=molecular_classification,eln_2017=eln_2017)


coeffs <- list()
dataTESTmodels <- list() 
bscmodels <- list()

bootstrap_iter <- 100

for (i in 1:length(prognosis_features)){    
    coeffs[names(prognosis_features)[i]] <- coef(cv.glmnet(as.matrix(dataTRAIN[,prognosis_features[[i]]]), as.matrix(yTRAIN), family="cox", alpha=0, nfolds=5, grouped=TRUE))
    RS <- as.matrix(dataTEST[,prognosis_features[[i]]]) %*% coeffs[[names(prognosis_features)[i]]]
    dataTESTmodels[names(prognosis_features)[i]] <- cbind(dataTEST$time, dataTEST$status, RS)
    bscmod <- NULL
    bscrsq <- NULL
    for(bs in 1:bootstrap_iter){
        set.seed(bs)  # so that all feature combination use same bootstrap data
        bssample <- sample(1:nrow(dataTEST), replace=TRUE)

        bscmod <- c(bscmod,Est.Cval(dataTESTmodels[[names(prognosis_features)[i]]][bssample,], tau=100)$Dhat)
        
    }
    bscmodels[[names(prognosis_features)[i]]] <- bscmod

}

write.table(bscmodels,"bootstrap_results_NEJM_cohort_eln_proposal.tsv")


################################## E. Subset of patients that received intensive treatment Comparison ELN vs Proposal ##################################

library(glmnet)
library(survC1)
library(survival)
options(warn=-1)


df_training <- read.table("aml_prognosis_updated.tsv")
df_training <- molecular_proposal_classification(df_training)

df_training$NEW_favorable <-  ifelse(df_training$molecular_classification=="NEW_favorable",1,0)
df_training$NEW_intermediate <-  ifelse(df_training$molecular_classification=="NEW_intermediate",1,0)
df_training$NEW_adverse <-  ifelse(df_training$molecular_classification=="NEW_adverse",1,0)
table(df_training$molecular_classification)

list_intensively_treated <- readRDS("list_intensively_treated.rds")
df_training_intense <- df_training[rownames(df_training) %in% list_intensively_treated,]
table(df_training_intense$molecular_classification)

df_training_intense$time <- df_training_intense$os
df_training_intense$status <- df_training_intense$os_status




########### Split the data into test and validation (75-25)
# For reproducibility, we use the same seed

IDTRAIN <- readRDS("list_NCRI_intense_Training.rds")
IDTEST <- readRDS("list_NCRI_intense_Testing.rds")

dataTRAIN <- df_training_intense[IDTRAIN,]
dataTEST <- df_training_intense[IDTEST,]




yTRAIN <- dataTRAIN[,c("time","status")]

molecular_classification <- c("NEW_favorable","NEW_intermediate","NEW_adverse")
eln_2017 <- c("eln_2017_favorable","eln_2017_intermediate","eln_2017_adverse")


prognosis_features<- list(molecular_classification=molecular_classification,eln_2017=eln_2017)


coeffs <- list()
dataTESTmodels <- list() 
bscmodels <- list()

bootstrap_iter <- 100

for (i in 1:length(prognosis_features)){
    coeffs[names(prognosis_features)[i]] <- coef(cv.glmnet(as.matrix(dataTRAIN[,prognosis_features[[i]]]), as.matrix(yTRAIN), family="cox", alpha=0, nfolds=5, grouped=TRUE))
    RS <- as.matrix(dataTEST[,prognosis_features[[i]]]) %*% coeffs[[names(prognosis_features)[i]]]
    dataTESTmodels[names(prognosis_features)[i]] <- cbind(dataTEST$time, dataTEST$status, RS)  
    bscmod <- NULL
    bscrsq <- NULL
    for(bs in 1:bootstrap_iter){
        set.seed(bs)  # so that all feature combination use same bootstrap data
        bssample <- sample(1:nrow(dataTEST), replace=TRUE)
        bscmod <- c(bscmod,Est.Cval(dataTESTmodels[[names(prognosis_features)[i]]][bssample,], tau=100)$Dhat)

    }
    bscmodels[[names(prognosis_features)[i]]] <- bscmod
}

write.table(bscmodels,"bootstrap_results_NCRI_cohort_eln_proposal_training_ONLY_intense.tsv")