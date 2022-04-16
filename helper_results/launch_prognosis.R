#.libPaths("/rigel/stats/users/yt2725/rpackages/")
library(glmnet)
library(parallel)
library(survival)
library(data.table)
library(CoxBoost)
library(randomForestSRC)
library(CoxHD)
library(survC1)
source('run_prognosis.R')

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

eln <- c(2,3,4)
age <- c(167)

all_gen <- c(5:88)
vect <- apply(X=df_final[,all_gen],2,FUN=function(x) 100*length(which(x==1))/dim(df_final)[1])
gen <- match(names(vect[vect>=2]),names(df_final))
              
all_cyto <- c(89:158)
vect <- apply(X=df_final[,all_cyto],2,FUN=function(x) 100*length(which(x==1))/dim(df_final)[1])
cyto <- match(names(vect[vect>=2]),names(df_final))
              
clin <- c(159:165)
demo <- c(166:167)

comp <-c(190:205)                        

clin_demo <-c(clin,demo)
              
all_gen_cyto <- c(all_gen,all_cyto)
              
all_gen_cyto_clin_demo <- c(all_gen_cyto,clin_demo)

comp_ITD <- c(comp,29)
              
comp_ITD_clin_demo <- c(comp_ITD,clin,demo)
              

y <- data.matrix(df_final[,c("os","os_status")])


prognosis_features<- list(all_gen_cyto=all_gen_cyto,clin_demo=clin_demo,all_gen_cyto_clin_demo=all_gen_cyto_clin_demo,comp_ITD_clin_demo=comp_ITD_clin_demo)
              
predictors <- c(rep(list(predictorGLM),1),rep(list(predictorRF),1),predictorBoost,predictorRFX)
str_predictors <-c(rep("CoxGLM",1),rep("RFS",1),"CoxBoost","RFX")
l_alpha <-seq(0,1,0.2)
l_ntree <- c(1050)
mc.cores <- 30
l_nodesize <- c(20)
for (i in 1:length(prognosis_features)){
    print("DONE")
    x <- data.matrix(df_final[,prognosis_features[[i]]])
    write.table(launch_prognosis(x=x,y=y,predictors=predictors,str_predictors=str_predictors,l_alpha=l_alpha,nrepeats=10,
                l_ntree=l_ntree,mc.cores=mc.cores,l_nodesize=l_nodesize),paste(names(prognosis_features)[i],".tsv",sep=""),quote=F,sep='\t')
    print("DONE")
    }