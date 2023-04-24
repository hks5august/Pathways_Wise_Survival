#Load required libraries
library(caret)
library(ggfortify)
library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(ggfortify)

#set working directory & seed 
setwd("/Users/kaurh8/Documents/CCGA_datasets/mRNA_693_samples/LGG/primary_vs_Recurrent/DGE/top_10pathways_survival/Apoptosis")
set.seed(7)


############################ Load Training data #######################
train_mat<-read.csv("../full_train_clin_23271.txt",header =TRUE, sep = "\t", check.names = FALSE)
head(train_mat[1:20],2)
dim(train_mat)
############################ Load validation data #######################
test1_mat<-read.csv("../full_test_clin_23271.txt",header =TRUE, sep = "\t", check.names = FALSE)
test2_mat<-read.csv("../full_ext_clin_23271.txt",header =TRUE, sep = "\t", check.names = FALSE)
dim(test1_mat)
dim(test2_mat)


######## Load file containing results (significant genes with p-value <0.05) from univariate analysis with beta coefficient  ##############################
all_features_results<-read.csv("../LGG_tr_Survival_results_for_genes.csv",header =TRUE, sep = "\t", check.names = FALSE)
dim(all_features_results)
head(all_features_results)

################################ Load selected features files (list of genes with ID as header) ##############################
#features<- read.csv(args[1], header =TRUE, sep = "\t", dec = ".",stringsAsFactors=FALSE, check.names = FALSE)
features<- read.csv("ApoptosisNew.tsv", header =TRUE, sep = "\t", dec = ".",stringsAsFactors=FALSE, check.names = FALSE)
dim(features)

######## Extract selected genes results from univariate analysis results, for the purpose of obtaining beta coefficient  ####################
sel_features_results<-all_features_results[row.names(all_features_results) %in% row.names(features),]

dim(sel_features_results)
head(sel_features_results)

#features1 <- features %>% dplyr::select(1,4,7:29)

sel_ftr_surv <- as.data.frame(sel_features_results[,1])
names(sel_ftr_surv ) <- c("ID")
head(sel_ftr_surv )
dim(sel_ftr_surv )

#### prepare training, test and external validation data with selected features having significant value in univariate analysis ##########
sel_train <- as.data.frame(train_mat[,colnames(train_mat) %in% c(sel_ftr_surv$ID), ])
head(sel_train,2)
sel_test <-as.data.frame(test1_mat[,colnames(test1_mat) %in% c(sel_ftr_surv$ID), ])
sel_ext_data <- as.data.frame(test2_mat[,colnames(test2_mat) %in% c(sel_ftr_surv$ID), ])

######### Make final files with the selected features (genes here)  and combine survival information ###########################################################
train_feature_mat<- cbind(train_mat["OS_month"],train_mat["Death_event"],sel_train)
head(train_feature_mat,2)
dim(train_feature_mat)
test1_feature_mat<- cbind(test1_mat["OS_month"],test1_mat["Death_event"],sel_test)
dim(test1_feature_mat)
test2_feature_mat<- cbind(test2_mat["OS_month"],test2_mat["Death_event"],sel_ext_data)
dim(test2_feature_mat)
#train_feature_mat<- cbind(train_mat["gender"],train_mat["age"],train_mat["OS"],train_mat["OS.time"],train_mat[,features])### Three blocks when files are loaded
#test1_feature_mat<- cbind(test_mat1["gender"],test_mat1["age"],test_mat1["OS"],test_mat1["OS.time"],test_mat1[,features])
#test2_feature_mat<- cbind(test_mat2["gender"],test_mat2["age"],test_mat2["OS"],test_mat2["OS.time"],test_mat2[,features])

################## results of selected features ### extract rows from the dataframe ####################
#sel_features_results<-all_features_results[row.names(all_features_results) %in% row.names(features),]

############ remove where OS.time=NA ############

train_feature_mat1<-subset(train_feature_mat,OS_month!="NA")
test1_feature_mat1<-subset(test1_feature_mat,OS_month!="NA")
test2_feature_mat1<-subset(test2_feature_mat,OS_month!="NA")

head(train_feature_mat1,2)
head(test1_feature_mat1,2)
head(test2_feature_mat1,2)
dim(train_feature_mat1)
dim(test1_feature_mat1)
dim(test2_feature_mat1)

#train_feature_mat2<- train_feature_mat1 %>% dplyr::select(1:2,10:23)
#test1_feature_mat2<- test1_feature_mat1 %>% dplyr::select(1:2,10:23)
#test2_feature_mat2<- test2_feature_mat1%>% dplyr::select(1:2,10:23)

#head(train_feature_mat2)

# save files with selected genes & survival information #########
write.table(cbind("ID"=rownames(train_feature_mat1), train_feature_mat1),file="sel_train.txt",sep="\t",quote=F, row.names=F)
write.table(cbind("ID"=rownames(test1_feature_mat1), test1_feature_mat1),file="sel_test1.txt",sep="\t",quote=F, row.names=F)
write.table(cbind("ID"=rownames(test2_feature_mat1), test2_feature_mat1),file="sel_test2.txt",sep="\t",quote=F, row.names=F)


#################### Training Data ##################################################################################

LUSC_C_tr=train_feature_mat1
tr <- LUSC_C_tr[3:ncol(LUSC_C_tr)]

E= length(tr)
E

PI_tr=0
for(i in seq(from=1, to=E ,by=1))
{
  PI_tr= PI_tr+(tr[,i]*sel_features_results[i,2])
}



LUSC_C_tr$PI<-PI_tr
head(LUSC_C_tr)

write.table(cbind("ID"=rownames(LUSC_C_tr), LUSC_C_tr),file="train_with_PI.txt",sep="\t",quote=F, row.names=F)
######################################## Survival Object ############################################

surv_object_tr <- Surv(time = LUSC_C_tr$OS_month, event = LUSC_C_tr$Death_event)

# ###survival analysis: fits cox ph model to find HR for PI
fit_tr <- survfit(surv_object_tr~(LUSC_C_tr$PI>mean(LUSC_C_tr$PI)), data=LUSC_C_tr)
summary(fit_tr)
fit.coxph_tr <- coxph(surv_object_tr ~(LUSC_C_tr$PI>mean(LUSC_C_tr$PI)), data=LUSC_C_tr)
summary(fit.coxph_tr)

tr_res <- summary(fit.coxph_tr)
coeff <- round(tr_res$coefficients[1],2)
HR <- round(tr_res$conf.int[1],2)
int1 <- round(tr_res$conf.int[3],2)
int2 <- round(tr_res$conf.int[4],2)
CI <- round (tr_res$concordance[1],2)

HR1 <- paste0(HR , " (",  int1,  " - ", int2, ") " )
HR1

tr_res <- cbind(coeff, HR1, CI)

#jpeg(file="test1_Beta_covar.jpeg", units="in", width=10, height=10, res=300)
jpeg(file="KM_plot_train_with_PI.jpeg", units="in", width=10, height=10, res=300)


ggsurvplot(fit_tr, data=LUSC_C_tr, pval=TRUE, 
           risk.table=TRUE, tables.height = 0.3, #add risk table & height
           xlab="Time in months", 
           risk.table.col="strata", break.time.by = 12,
           conf.int = F, censor = TRUE,
           surv.median.line = "hv", # Add medians survival
           #palette = c("red","blue"),#add desired color
           size=1.5, font.tickslab = c(12,  "black"), font.y = c(14,  "black"), 
           font.x = c(14, "black"), font.legend = c(14, "bold", "black"), font.label = c(14, "bold", "black"))

dev.off()


################################# Test1 #############################################################
dim(LUSC_C_te1)
dim(sel_features_results)
LUSC_C_te1=test1_feature_mat1

te1 <- LUSC_C_te1[3:ncol(LUSC_C_te1)]

E= length(te1)
E

PI_te1=0
for(i in seq(from=1, to=E ,by=1))
{
  PI_te1= PI_te1+(te1[,i]*sel_features_results[i,2])
}
PI_te1


LUSC_C_te1$PI <- PI_te1
dim(LUSC_C_te1)
head(LUSC_C_te1)

write.table(cbind("ID"=rownames(LUSC_C_te1), LUSC_C_te1),file="test1_with_PI.txt",sep="\t",quote=F, row.names=F)
######################################## Survival Object ############################################

surv_object_te1 <- Surv(time = LUSC_C_te1$OS_month, event = LUSC_C_te1$Death_event)

# ###survival analysis: fits cox ph model to find HR for PI
fit_te1 <- survfit(surv_object_te1~(LUSC_C_te1$PI>mean(LUSC_C_te1$PI)), data=LUSC_C_te1)
summary(fit_te1)
fit.coxph_te1 <- coxph(surv_object_te1 ~(LUSC_C_te1$PI>mean(LUSC_C_te1$PI)), data=LUSC_C_te1)
#ggsurvplot(fit)
# write.csv(LUSC_C, file = "/Users/sherrybhalla/Desktop/LUNG-Cancer/LUNG-Apop-HRmean-PImean.csv",row.names = FALSE)

te1_res <- summary(fit.coxph_te1)
coeff <- round(te1_res$coefficients[1],2)
HR <- round(te1_res$conf.int[1],2)
int1 <- round(te1_res$conf.int[3],2)
int2 <- round(te1_res$conf.int[4],2)
CI <- round (te1_res$concordance[1],2)

HR1 <- paste0(HR , " (",  int1,  " - ", int2, ") " )
HR1

te1_res1 <- cbind(coeff, HR1, CI)

#jpeg(file="test1_Beta_covar.jpeg", units="in", width=10, height=10, res=300)
jpeg(file="KM_plot_test1_with_PI.jpeg", units="in", width=10, height=10, res=300)


ggsurvplot(fit_te1, data=LUSC_C_te1, pval=TRUE, 
           risk.table=TRUE, tables.height = 0.3, #add risk table & height
           xlab="Time in months", 
           risk.table.col="strata", break.time.by = 12,
           conf.int = F, censor = TRUE,
           surv.median.line = "hv", # Add medians survival
           #palette = c("red","blue"),#add desired color
           size=1.5, font.tickslab = c(12,  "black"), font.y = c(14,  "black"), 
           font.x = c(14, "black"), font.legend = c(14, "bold", "black"), font.label = c(14, "bold", "black"))

dev.off()

################################# Test2 #############################################################

LUSC_C_te2=test2_feature_mat1


te2 <- LUSC_C_te2[3:ncol(LUSC_C_te2)]

E= length(te2)
E



PI_te2=0
for(i in seq(from=3, to=E ,by=1))
{
  PI_te2= PI_te2+(te2[,i]*sel_features_results[i,2])
}

LUSC_C_te2$PI<-PI_te2

write.table(cbind("ID"=rownames(LUSC_C_te2), LUSC_C_te2),file="test2_with_PI.txt",sep="\t",quote=F, row.names=F)
######################################## Survival Object ############################################

surv_object_te2 <- Surv(time = LUSC_C_te2$OS_month, event = LUSC_C_te2$Death_event)

# ###survival analysis: fits cox ph model to find HR for PI
fit_te2 <- survfit(surv_object_te2~(LUSC_C_te2$PI>mean(LUSC_C_te2$PI)), data=LUSC_C_te2)
summary(fit_te2)
fit.coxph_te2 <- coxph(surv_object_te2 ~(LUSC_C_te2$PI>mean(LUSC_C_te2$PI)), data=LUSC_C_te2)
#ggsurvplot(fit)
# write.csv(LUSC_C, file = "/Users/sherrybhalla/Desktop/LUNG-Cancer/LUNG-Apop-HRmean-PImean.csv",row.names = FALSE)

te2_res <- summary(fit.coxph_te2)
coeff <- round(te2_res$coefficients[1],2)
HR <- round(te2_res$conf.int[1],2)
int1 <- round(te2_res$conf.int[3],2)
int2 <- round(te2_res$conf.int[4],2)
CI <- round (tr_res$concordance[1],2)

HR1 <- paste0(HR , " (",  int1,  " - ", int2, ") " )
HR1

te2_res1 <- cbind(coeff, HR1, CI)

#jpeg(file="test1_Beta_covar.jpeg", units="in", width=10, height=10, res=300)
jpeg(file="KM_plot_test2_with_PI.jpeg", units="in", width=10, height=10, res=300)


ggsurvplot(fit_te2, data=LUSC_C_te2, pval=TRUE, 
           risk.table=TRUE, tables.height = 0.3, #add risk table & height
           xlab="Time in months", 
           risk.table.col="strata", break.time.by = 12,
           conf.int = F, censor = TRUE,
           surv.median.line = "hv", # Add medians survival
           #palette = c("red","blue"),#add desired color
           size=1.5, font.tickslab = c(12,  "black"), font.y = c(14,  "black"), 
           font.x = c(14, "black"), font.legend = c(14, "bold", "black"), font.label = c(14, "bold", "black"))

dev.off()


beta_cov <- rbind(tr_res1, te1_res1,  te2_res1)
beta_cov

write.table(beta_cov , file = "Pathway_PI_results.txt",row.names = T, col.names=T, sep="\t")


###################################################################################
