################################## A. R squared training and validation whole cohort ##################################

# Bootstrap SG and NCRI 100 times for R squared
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

library(parallel)
library(survival)
library(ggpubr)
df_training <- read.table("aml_prognosis_updated.tsv")
df_training_proposal <- molecular_proposal_classification(df_training)

df_validation <- read.table("full_data_validation.tsv")
df_validation_proposal <- molecular_proposal_classification(df_validation)

table(df_training_proposal$molecular_classification)
table(df_validation_proposal$molecular_classification)

set.seed(17)
res_bootstrap <- data.frame('rsq' = numeric(),
                  'concordance' = numeric())
n_exp = 100
folds_training <- list()
folds_validation <- list()

for (i in seq(n_exp)) {
    folds_training[[i]] <- sample(1:nrow(df_training_proposal), 0.75 * nrow(df_training_proposal), replace = TRUE)
    folds_validation[[i]] <- sample(1:nrow(df_validation_proposal), 0.75 * nrow(df_validation_proposal), replace = TRUE)
}
nexp = length(folds_training)
print("Start Bootstrapping")
rescv = mclapply(seq(nexp),
               FUN=function(iexp) {
                   set.seed(17)
                   x_training <- df_training_proposal[folds_training[[iexp]],]
                   x_validation <- df_validation_proposal[folds_validation[[iexp]],]
                   tmp <- c(training_rsquared_proposal=summary(coxph(Surv(os,os_status)~molecular_classification,x_training))$rsq[[1]],
                            training_rsquared_eln=summary(coxph(Surv(os,os_status)~eln_2017,x_training))$rsq[[1]],
                            validation_rsquared_proposal=summary(coxph(Surv(OS,OS_Status)~molecular_classification,x_validation))$rsq[[1]],
                            validation_rsquared_eln=summary(coxph(Surv(OS,OS_Status)~eln_2017,x_validation))$rsq[[1]])
               },
               mc.cores=16
               )

df_Rsquared <- data.frame("Dataset"=c(rep("Training Cohort",200),rep("Validation Cohort",200)),"Model"=c(rep("Model",100),rep("ELN",100),rep("Model",100),rep("ELN",100))
                  ,"R_squared"=rep(0,400))
data_comparison <- data.frame(matrix(unlist(rescv), nrow=length(rescv), byrow=T))
colnames(data_comparison) <- names(rescv[[1]])
df_Rsquared[1:100,"R_squared"] <- data_comparison$training_rsquared_proposal
df_Rsquared[101:200,"R_squared"] <- data_comparison$training_rsquared_eln

df_Rsquared[201:300,"R_squared"] <- data_comparison$validation_rsquared_proposal
df_Rsquared[301:400,"R_squared"] <- data_comparison$validation_rsquared_eln

anno_df_rsq <- compare_means(R_squared~Model, group.by = "Dataset", data = df_Rsquared)

write.table(df_Rsquared[df_Rsquared$Dataset=="Training Cohort",],"results_NCRI_cohort_eln_proposal_RSQUARED.tsv")
write.table(df_Rsquared[df_Rsquared$Dataset=="Validation Cohort",],"results_NEJM_cohort_eln_proposal_RSQUARED.tsv")




################################## B. R squared training and validation whole cohort ##################################

# Bootstrap SG and NCRI 100 times for R squared

library(parallel)
library(survival)
library(ggpubr)
df_training <- read.table("aml_prognosis_updated.tsv")
df_training <- molecular_proposal_classification(df_training)
table(df_training$molecular_classification)


list_intensively_treated <- readRDS("list_intensively_treated.rds")
df_training_intense <- df_training[rownames(df_training) %in% list_intensively_treated,]


set.seed(17)
res_bootstrap <- data.frame('rsq' = numeric(),
                  'concordance' = numeric())
n_exp = 100
folds_training <- list()

for (i in seq(n_exp)) {
    folds_training[[i]] <- sample(1:nrow(df_training_intense ), 0.75 * nrow(df_training_intense ), replace = TRUE)
}
nexp = length(folds_training)
print("Start Bootstrapping")
rescv = mclapply(seq(nexp),
               FUN=function(iexp) {
                   set.seed(17)
                   x_training <- df_training_intense [folds_training[[iexp]],]
                   tmp <- c(training_rsquared_proposal=summary(coxph(Surv(os,os_status)~molecular_classification,x_training))$rsq[[1]],
                            training_rsquared_eln=summary(coxph(Surv(os,os_status)~eln_2017,x_training))$rsq[[1]])
               },
               mc.cores=16
               )

df_Rsquared <- data.frame("Dataset"=rep("Training Cohort",200),"Model"=c(rep("Model",100),rep("ELN",100)),"R_squared"=rep(0,200))
data_comparison <- data.frame(matrix(unlist(rescv), nrow=length(rescv), byrow=T))
colnames(data_comparison) <- names(rescv[[1]])
df_Rsquared[1:100,"R_squared"] <- data_comparison$training_rsquared_proposal
df_Rsquared[101:200,"R_squared"] <- data_comparison$training_rsquared_eln

anno_df_rsq <- compare_means(R_squared~Model, group.by = "Dataset", data = df_Rsquared)

write.table(df_Rsquared,"results_NCRI_cohort_eln_proposal_RSQUARED_ONLY_INTENSE.tsv")