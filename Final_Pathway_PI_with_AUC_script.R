#Load required libraries
library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(ggfortify)
library(survivalROC)
library(MASS)
args <- commandArgs(TRUE)

set.seed(7)

path <- paste0("/Users/kaurh8/Documents/Pathways_Survival_analysis_all_cancer/")

#set working directory
setwd(path)
getwd()
cancer_list <- read.table("list2",header =T, sep = "\t",  row.names=1,  check.names = FALSE)
folder <- as.character(row.names(cancer_list))
folder




for (c in 1:length(folder))
{
  setwd(paste0(path,folder[c], "/"))
  getwd()


#setwd("/Users/kaurh8/Documents/Pathways_Survival_analysis_all_cancer/LIHC/")
#}
############################ Load Training data #######################
#data<-read.csv("LIHC.csv",header =TRUE, sep = ",", row.names=1,  check.names = FALSE)
data<-read.table(paste0(path,folder[c],".csv"),header =TRUE, sep = ",", row.names=1,  check.names = FALSE)
dim(data)

#extract clinical data
Clin_data <- data[1:14]
dim(Clin_data)
Clin_data1 <- Clin_data %>% mutate(OS_month = round(OS.time/30.417, digit=0))

#extract expression data
exp_data <- data[15:ncol(data)]
dim(exp_data)

#combine new clinical data with exp data
train_mat <- cbind(Clin_data1, exp_data)


######## Load file containing results (significant genes with p-value <0.05) from univariate analysis with beta coefficient  ##############################
all_features_results<-read.csv("Significant_Survival_results_for_genes.txt",header =TRUE, sep = "\t", check.names = FALSE)
dim(all_features_results)
head(all_features_results)

pathway_list <- read.table("../3pathways_list",header =T, sep = "\t",  row.names=1,  check.names = FALSE)

#path <- "/Users/kaurh8/Documents/Pathways_Survival_analysis_all_cancer/LIHC/"
subfolder <- as.character(row.names(pathway_list))
subfolder

#create directory for each pathway within each cancer then create prognostic index analysis for each pathway 

#create file containing results for all pathways
write.table(cbind("ID","Beta","HR_with_CI","Concordance","P-value"),
            file=paste0(path,folder[c], "/", "Pathways_PI_Survival_results.txt"),row.names=F,col.names=F,sep = '\t');

for (j in 1:length(subfolder))
  {

  folder1<-dir.create(paste0(path,folder[c], "/",subfolder[j]))

a<- paste0(folder1, "/")
print(a)  
  
  setwd(paste0(path,folder[c], "/",subfolder[j]))
getwd()
  
  
  path2 <- "/Users/kaurh8/Documents/Pathways_Survival_analysis_all_cancer/All_pathways/"
  
#load pathway file containing list of genes
 #features<- read.csv("ApoptosisNew.tsv", header =TRUE, sep = "\t", dec = ".",stringsAsFactors=FALSE, check.names = FALSE)
  #features<- read.csv(args[1], header =TRUE, sep = "\t", dec = ".",stringsAsFactors=FALSE, check.names = FALSE)
  features<- read.csv(paste0(path2,subfolder[j], ".txt"), header =TRUE, sep = "\t", dec = ".",stringsAsFactors=FALSE, check.names = FALSE)
  dim(features)
#save pathway file within directory  
write.table(cbind("ID"=rownames(features ), features ),file="features.txt",sep="\t",quote=F, row.names=F)
  
  ######## Extract selected genes results from univariate analysis results, for the purpose of obtaining beta coefficient  ####################
  sel_features_results<-all_features_results[row.names(all_features_results) %in% row.names(features),]
  
  dim(sel_features_results)
  head(sel_features_results)
  
  
  sel_ftr_surv <- as.data.frame(sel_features_results[,1])
  names(sel_ftr_surv ) <- c("ID")
  head(sel_ftr_surv )
  dim(sel_ftr_surv )
  
  #### prepare training, test and external validation data with selected features having significant value in univariate analysis ##########
  sel_train <- as.data.frame(train_mat[,colnames(train_mat) %in% c(sel_ftr_surv$ID), ])
  head(sel_train,2)
  
  ######### Make final files with the selected features (genes here)  and combine survival information ###########################################################
  train_feature_mat<- cbind(train_mat["OS_month"],train_mat["OS"],sel_train)
  head(train_feature_mat,2)
  dim(train_feature_mat)
  
  ############ remove where OS.time=NA ############
  train_feature_mat1<-subset(train_feature_mat,OS_month!="NA")
  
  
  # save files with selected genes & survival information #########
  write.table(cbind("ID"=rownames(train_feature_mat1), train_feature_mat1),file="sel_train.txt",sep="\t",quote=F, row.names=F)
  
#Create prognostic index for pathway 
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
  
#save selected data with PI value  
write.table(cbind("ID"=rownames(LUSC_C_tr), LUSC_C_tr),file="train_with_PI.txt",sep="\t",quote=F, row.names=F)
  
  tr_PI <- as.data.frame(LUSC_C_tr$PI)
  rownames(tr_PI) <- rownames(LUSC_C_tr)
  colnames(tr_PI) <- c("PI")
  write.table(cbind("ID"=rownames(tr_PI ), tr_PI ),file="tr_PI",sep="\t",quote=F, row.names=F)
  
  
  ######################################## Survival Object ############################################
  surv_object_tr <- Surv(time = LUSC_C_tr$OS_month, event = LUSC_C_tr$OS)
  dim(surv_object_tr )
  dim(LUSC_C_tr$PI)
  
  dim(tr)
  head(tr_PI)
  mean(tr_PI$PI)
  
  # ###survival analysis: fits cox ph model to find HR for PI
  fit_tr <- survfit(surv_object_tr~(tr_PI$PI>mean(tr_PI$PI)), data=tr_PI)
  #fit_tr <- survfit(surv_object_tr~(LUSC_C_tr$PI>mean(LUSC_C_tr$PI)), data=LUSC_C_tr)
  summary(fit_tr)

  #fit.coxph_tr <- coxph(surv_object_tr ~(LUSC_C_tr$PI>mean(LUSC_C_tr$PI)), data=LUSC_C_tr)
  fit.coxph_tr <- coxph(surv_object_tr ~(tr_PI$PI>mean(tr_PI$PI)), data=tr_PI)
  summary(fit.coxph_tr)
  
  
  tr_res <- summary(fit.coxph_tr)
  coeff <- round(tr_res$coefficients[1],2)
  HR <- round(tr_res$conf.int[1],2)
  int1 <- round(tr_res$conf.int[3],2)
  int2 <- round(tr_res$conf.int[4],2)
  CI <- round (tr_res$concordance[1],2)
  pval <- tr_res$sctest[3]
  
  HR1 <- paste0(HR , " (",  int1,  " - ", int2, ") " )
  HR1
  
  tr_res1 <- cbind(coeff, HR1, CI, pval)
  tr_res1 
#save results as a file  
write.table(cbind("ID"=rownames(tr_res1), tr_res1),file="tr_res1.txt",sep="\t",quote=F, row.names=F)

#KM plot
#set theme 
centr = theme_survminer() + theme(plot.title = element_text(hjust = 0.5, face = "bold"))

#save KM survival plot
#jpeg(file="KM_plot.jpeg", units="in", width=15, height=10, res=300)
pp <-  ggsurvplot(fit_tr, data=tr_PI, pval=TRUE,
                  risk.table=TRUE, tables.height = 0.3, #add risk table & height
                  xlab="Time in Months",
                  risk.table.col="strata", break.time.by = 12, fontsize = 3,
                  conf.int = F, censor = TRUE,
                  title= paste0("KM plot based on PI of ", pathway , " in ", cancer),
                  ggtheme=centr,
                   surv.median.line = "hv", # Add medians survival
                  #palette = c("red","blue"),#add desired color
                  size=1.2, font.tickslab = c(10,  "black"), font.y = c(10,  "black"),
                  font.x = c(10, "black"), font.legend = c(10, "bold", "black"), font.label = c(10, "bold", "black"))
#pp
jpeg("KM_plot.jpg", units="in", width=15, height=10, res=300)
print(pp, newpage = FALSE)
dev.off()


#Create ROC plot for 1-,3- and 5-years survival time prediction
  
  #Calculated  Lmbda min value = 0.002102803
  
  tr_roc1 <- survivalROC(Stime        = LUSC_C_tr$OS_month,
                         status       = LUSC_C_tr$OS,
                         marker       = tr_PI$PI,
                         predict.time = 12,
                         method       = "NNE", lambda = 0.0021, span = NULL, window ="symmetric")
  
  
  tr_roc3 <- survivalROC(Stime        = LUSC_C_tr$OS_month,
                         status       = LUSC_C_tr$OS,
                         marker       = tr_PI$PI,
                         predict.time = 36,
                         method       = "NNE", lambda = 0.0021, span = NULL, window ="symmetric")
  
  tr_roc5 <- survivalROC(Stime        = LUSC_C_tr$OS_month,
                         status       = LUSC_C_tr$OS,
                         marker       = tr_PI$PI,
                         predict.time = 60,
                         method       = "NNE", lambda = 0.0021, span = NULL, window ="symmetric")
  
  jpeg(file="ROC.jpeg", units="in", width=10, height=10, res=300)
  plot(tr_roc1$FP, tr_roc1$TP, type="l", xlim=c(0,1), ylim=c(0,1),
       xlab="FP", col="red", 
       ylab="TP",main= "AUC Curve for Survival Prediction")
  lines(tr_roc3$FP, tr_roc3$TP, type="l", lty=2, col="blue")
  lines(tr_roc5$FP, tr_roc5$TP, type="l", lty=2, col="green")
  legend(0.5,0.5, legend=c( paste( "1 Year AUC = ",round(tr_roc1$AUC,2)),  paste("3 Years AUC = ",round(tr_roc3$AUC,2)),  paste("5 Years AUC = ",round(tr_roc5$AUC,2))), col =c ("red","blue", "green"), lty=c(1,2), bty="n")
  
  abline(0,1)
  
  dev.off()
 
  write.table(cbind(paste0(subfolder[j]),coeff, HR1, CI, pval),
              file=paste0(path,folder[c], "/", "Pathways_PI_Survival_results.txt"),row.names=F,col.names=F,sep = '\t',append = T);#output file
  
  
}


}
