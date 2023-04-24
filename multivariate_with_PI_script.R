
#Load required libraries
library(ggfortify)
library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(ggfortify)

#library(pca3d)


#setwd("/Users/kaurh8/Documents/CCGA_datasets/mRNA_693_samples/LGG/LGG/survival")
#setwd("/Users/kaurh8/Documents/CCGA_datasets/mRNA_693_samples/LGG/primary_vs_Recurrent/DGE/top_10pathways_survival/Apoptosis")
set.seed(7)



tr_clin_new <- read.table("tr_clin_PI", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
te1_clin_new <- read.table("te1_clin_PI", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
te2_clin_new <- read.table("te2_clin_PI", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

dim(tr_clin_new)
dim(te1_clin_new)
dim(te2_clin_new)




# create survival object
surv_object_p_tr <- Surv(time = as.numeric(tr_clin_new$OS_month), event = tr_clin_new$Death_event)
#Cox multivariate model
cox_multivariate_tr <- coxph(surv_object_p_tr ~  Gender + Class + Grade + as.numeric(Age) + IDH_mutation_status + codel_1p19q_status + MGMTp_meth_status + RT + Chemo_TMZ  + PI,  data=tr_clin_new  )
summary(cox_multivariate_tr )
jpeg(file="train_multivariate.jpeg", units="in", width=10, height=10, res=300)
#ggforest(cox_multivariate,data=LGG_clin )
ggforest(cox_multivariate_tr,data=tr_clin_new)
dev.off()



# create survival object
surv_object_p_te1 <- Surv(time = as.numeric(te1_clin_new$OS_month), event = te1_clin_new$Death_event)
#Cox multivariate model
cox_multivariate_te1 <- coxph(surv_object_p_te1 ~  Gender + Class + Grade + as.numeric(Age) + IDH_mutation_status + codel_1p19q_status + MGMTp_meth_status + RT + Chemo_TMZ  + PI,  data=te1_clin_new  )
summary(cox_multivariate_te1 )
jpeg(file="Test1_multivariate.jpeg", units="in", width=10, height=10, res=300)
#ggforest(cox_multivariate,data=LGG_clin )
ggforest(cox_multivariate_te1,data=te1_clin_new)
dev.off()



# create survival object
surv_object_p_te2 <- Surv(time = as.numeric(te2_clin_new$OS_month), event = te2_clin_new$Death_event)
#Cox multivariate model
cox_multivariate_te2 <- coxph(surv_object_p_te2 ~  Gender + Class + Grade + as.numeric(Age) + IDH_mutation_status + codel_1p19q_status + MGMTp_meth_status + RT + Chemo_TMZ  + PI,  data=te2_clin_new  )
summary(cox_multivariate_te2 )
jpeg(file="Test2_multivariate.jpeg", units="in", width=10, height=10, res=300)
#ggforest(cox_multivariate,data=LGG_clin )
ggforest(cox_multivariate_te2,data=te2_clin_new)
dev.off()

