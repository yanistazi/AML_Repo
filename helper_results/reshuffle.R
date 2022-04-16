## Launch importance to compare the important features coming out from comp genes cytos on Columbia's cluster
#.libPaths("/rigel/stats/users/yt2725/rpackages/")
library(glmnet)
library(parallel)
library(survival)
library(data.table)
library(CoxBoost)
library(randomForestSRC)
library(CoxHD)
options(warn=-1)
source('run_prognosis.R')
source("feature_importance.R")



### Features that we can use
###-----------------------------------------------------------------------------
df_final <- read.table("aml_prognosis_updated.tsv")

normalized<-function(y) {

  x<-y[!is.na(y)]

  x<-(x - min(x)) / (max(x) - min(x))

  y[!is.na(y)]<-x

  return(y)
}

df_final$gender <- factor(df_final$gender,labels=c(0,1))
df_final$ahd <- factor(df_final$ahd,labels=c(0,1))
df_final$perf_status <- factor(df_final$perf_status,labels=c(0,1,2,3,4))
df_final$secondary <- factor(df_final$secondary,labels=c(1,2,3))
df_final[,c("age","bm_blasts","wbc","hb","plt")] <- apply(df_final[,c("age","bm_blasts","wbc","hb","plt")],2,normalized) 

all_gen <- c(5:88)
vect <- apply(X=df_final[,all_gen],2,FUN=function(x) 100*length(which(x==1))/dim(df_final)[1])
gen <- match(names(vect[vect>=2]),names(df_final))
              
all_cyto <- c(89:158)
vect <- apply(X=df_final[,all_cyto],2,FUN=function(x) 100*length(which(x==1))/dim(df_final)[1])
cyto <- match(names(vect[vect>=2]),names(df_final))
              
clin <- c(159:165)
demo <- c(166:167)
demo_without_age <-c(166)
comp <- c(190:205)
eln <- c(2,3,4)            

all_gen_cyto <- c(all_gen,all_cyto)
                          
comp_gen_cyto <- c(comp,gen,cyto)
comp_all_gen_cyto <- c(comp,all_gen_cyto)
eln_all_gen_cyto <- c(eln,all_gen_cyto)
eln_all_gen <- c(eln,all_gen)
eln_all_gen_cyto_clin_demo <- c(eln,all_gen_cyto,clin,demo)

y <- data.matrix(df_final[,c("os","os_status")])



prognosis_features<- list(all_gen=all_gen,comp_all_gen_cyto=comp_all_gen_cyto,comp_gen_cyto=comp_gen_cyto,comp=comp,eln_all_gen_cyto=eln_all_gen_cyto,
                          eln_all_gen=eln_all_gen,eln_all_gen_cyto_clin_demo=eln_all_gen_cyto_clin_demo)          
##---------------------------------------------------------------------------------PREPARING MODELS and ALGOS
                         

nrepeats=5
seed=1234
mc.cores=30
npermutations=4
nfolds=5

algorithms<-c(algo_Lasso, algo_Ridge, algo_Elastic_net,  algo_RFX, algo_RFS, algo_Cox)
predictors<-c(predictor_Lasso, predictor_Ridge, predictor_Elastic_net,  predictor_RFX, predictor_RFS,  predictor_Cox)
algo_names<-c('Lasso','Ridge','Elastic_net','RFX','RFS','Cox')

response <- data.matrix(df_final[,c("os","os_status")])
colnames(response) <- c("time","status")


#---------------------------------------------------------------------------------PREPARING MODELS and ALGOS

for (j in 1:length(prognosis_features)){
    print(names(prognosis_features[j]))
    res_data <- data.frame('feature'=character(),'ref_CI'=numeric(),'permuted_CI'=numeric(),'algo'=character(),'model'=character())
    for(i in 1:length(algorithms)){
        design <- data.matrix(data.frame(df_final[,prognosis_features[[j]]]))      
        tmp <- runCV_CI_with_test(response=response, design=design,
              nfolds=nfolds, nrepeats=nrepeats, seed=seed, mc.cores=mc.cores, features=colnames(design), npermutations=npermutations, 
                                  algorithm=algorithms[i][[1]], predictor=predictors[i][[1]])
        tmp$algo<-algo_names[i]
        tmp$model <- names(prognosis_features[j])
        res_data <- rbind(res_data,tmp)
    }
    write.table(res_data,paste(names(prognosis_features)[j],".tsv",sep="_reshuffle_importance"),quote=F,sep='\t')
}
              
                           
              
