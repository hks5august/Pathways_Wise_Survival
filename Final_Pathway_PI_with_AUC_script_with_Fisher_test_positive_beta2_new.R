#usage: nohup Rscript Final_Pathway_PI_with_AUC_script_copy.R cancer_list pathways_list &
#requirement   folder for each cancer containing two files  (files gene expression, where samples are in row and genes in column; univariate survival analysis results, where for each gene we have beta coefficient )
#list of cancer and list of pathways; folder containing files with all genes of  each pathway
#Load required libraries
library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(ggfortify)
library(survivalROC)
library(MASS)
library(glmnet)

args <- commandArgs(TRUE)

set.seed(7)

#path wher have input data all files and folders
#path <- paste0("/Users/kaurh8/Documents/Pathways_Survival_analysis_all_cancer/")
path <- paste0("/data/kaurh8/Pathways_PI/Final_Pathway_Folder/")
#set working directory
setwd(path)
getwd()

#provide list of cancer type
cancer_list <- read.table(paste0(path, "cancer_list4"),header =T, sep = "\t",  row.names=1,  check.names = FALSE)
#cancer_list <- read.table(paste0(path, "Final_cancer_list"),header =T, sep = "\t",  row.names=1,  check.names = FALSE)
#cancer_list <- read.table(args[1],header =T, sep = "\t",  row.names=1,  check.names = FALSE)
#cancer_list <- read.table(paste0(path, "cancer_list"),header =T, sep = "\t",  row.names=1,  check.names = FALSE)
folder <- as.character(row.names(cancer_list))
folder
#print("folder:" folder)


for (c in 1:length(folder))
{
  setwd(paste0(path,folder[c], "/"))
 getwd()

path_c <- paste0(path,folder[c], "/")
#print("cancer folder", path_c)
path_c
setwd(path_c)
getwd()
#setwd("/Users/kaurh8/Documents/Pathways_Survival_analysis_all_cancer/LIHC/")
#}
############################ Load Training data #######################
#data<-read.csv("LIHC.csv",header =TRUE, sep = ",", row.names=1,  check.names = FALSE)
#data<-read.table(paste0(path,folder[c],".csv"),header =TRUE, sep = "\t", row.names=1,  check.names = FALSE)
data<-read.table(paste0(path_c,folder[c],".csv"),header =TRUE, sep = "\t", row.names=1,  check.names = FALSE)
dim(data)

#extract clinical data
#Clin_data <- data[1:14]
Clin_data <- data[1:15]
dim(Clin_data)
#Days into months
#Clin_data1 <- Clin_data %>% mutate(OS_months = round(OS.time/30.417, digit=0))
#Years into months
Clin_data1 <- Clin_data %>% mutate(OS_months = round(OS.time * 12, digit=0))
#extract expression data
exp_data <- data[16:ncol(data)]
dim(exp_data)

#combine new clinical data with exp data
#train_mat <- cbind(Clin_data1, exp_data)
data2 <- cbind(Clin_data1, exp_data)

#remove samples with NA values of OS time 
train_mat <- subset(data2, OS.time!="NA")
train_mat <- subset(train_mat, OS.time >0)

######## Load file containing results (significant genes with p-value <0.05) from univariate analysis with beta coefficient  ##############################
all_features_results<-read.table(paste0(path_c,"Significant_Survival_results_for_genes.txt"),header =TRUE, sep = "\t", check.names = FALSE)
#all_features_results <- read.table(paste0(path_c,"DGE_Prog_Significant_Survival_results_for_genes.txt"), header =TRUE, sep = "\t", check.names = FALSE)

pos_features <- subset(all_features_results, Beta > 0 )
dim(all_features_results)
head(all_features_results)

#provide pathway list
#pathway_list <- read.table("../3pathways_list",header =T, sep = "\t",  row.names=1,  check.names = FALSE)
#pathway_list <- read.table(args[2],header =T, sep = "\t",  row.names=1,  check.names = FALSE)
#pathway_list <- read.table(paste0(path, "pathways_list"),header =T, sep = "\t",  row.names=1,  check.names = FALSE)
pathway_list <- read.table(paste0(path, "final_pathways_list"),header =T, sep = "\t",  row.names=1,  check.names = FALSE)
#pathway_list <- read.table(paste0(path, "final_pathways_list "),header =T, sep = "\t",  row.names=1,  check.names = FALSE)

#path <- "/Users/kaurh8/Documents/Pathways_Survival_analysis_all_cancer/LIHC/"
subfolder <- as.character(row.names(pathway_list))
subfolder

#create directory for each pathway within each cancer then create prognostic index analysis for each pathway 

#create file containing results for all pathways
write.table(cbind("Pathway","Beta","HR_with_CI","Concordance","P-value", "Total_path_genes" , "Path_sig_Yes", "Percentage", "p.value", "odd_ratio"),
            file=paste0(path,folder[c], "/", "Pos_Pathways_PI_Survival_results.txt"),row.names=F,col.names=F,sep = '\t');

for (j in 1:length(subfolder))
  {
tryCatch({

  folder1<-dir.create(paste0(path,folder[c], "/","Pos_", subfolder[j]))

a<- paste0(folder1, "/")
print(a)  
  
  setwd(paste0(path,folder[c], "/", "Pos_", subfolder[j]))
getwd()
  
  
 # path2 <- "/Users/kaurh8/Documents/Pathways_Survival_analysis_all_cancer/All_pathways/"
 path2 <- "/data/kaurh8/Pathways_PI/Final_Pathway_Folder/All_pathways/" 
#load pathway file containing list of genes
 #features<- read.csv("ApoptosisNew.tsv", header =TRUE, sep = "\t", dec = ".",stringsAsFactors=FALSE, check.names = FALSE)
  #features<- read.csv(args[1], header =TRUE, sep = "\t", dec = ".",stringsAsFactors=FALSE, check.names = FALSE)
  features<- read.csv(paste0(path2,subfolder[j], ".txt"), header =TRUE, sep = "\t", dec = ".",stringsAsFactors=FALSE, check.names = FALSE)
  dim(features)
#save pathway file within directory  
write.table(cbind("ID"=rownames(features ), features ),file="features.txt",sep="\t",quote=F, row.names=F)
  
  ######## Extract selected genes results from univariate analysis results, for the purpose of obtaining beta coefficient  ####################
 # sel_features_results<-all_features_results[row.names(all_features_results) %in% row.names(features),]
 
  sel_features_results<-pos_features[row.names(pos_features) %in% row.names(features),] 
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
  train_feature_mat<- cbind(train_mat["OS.time"],train_mat["OS"],sel_train)
  head(train_feature_mat,2)
  dim(train_feature_mat)
  
  ############ remove where OS.time=NA ############
  train_feature_mat1<-subset(train_feature_mat,OS.time!="NA")
  train_feature_mat1<-subset(train_feature_mat1,OS.time!=0)
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
    PI_tr= PI_tr+((tr[,i])*(sel_features_results[i,2]))
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
  surv_object_tr <- Surv(time = LUSC_C_tr$OS.time, event = LUSC_C_tr$OS)
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
 
#save KM survival plot 
  #jpeg(file= paste0(path, subfolder[j], "/", "KM_plot_train_with_PI.jpeg"), units="in", width=10, height=10, res=300)
  #paste0(path, subfolder[j], "/", "KM_plot_train_with_PI.jpeg")
#jpeg(file="KM_plot.jpeg", units="in", width=10, height=10, res=300)
pp <-  ggsurvplot(fit_tr, data=tr_PI, pval=TRUE,
             risk.table=TRUE, tables.height = 0.3, #add risk table & height
             xlab="Time in Days",
             risk.table.col="strata", break.time.by = 12,
             conf.int = F, censor = TRUE,
             surv.median.line = "hv", # Add medians survival
             #palette = c("red","blue"),#add desired color
             size=1.5, font.tickslab = c(12,  "black"), font.y = c(14,  "black"),
             font.x = c(14, "black"), font.legend = c(14, "bold", "black"), font.label = c(14, "bold", "black"))
#pp
jpeg("KM_plot.jpg", units="in", width=10, height=10, res=300)
print(pp, newpage = FALSE)
dev.off()

###################### Fisher Test ###########
# calculate number of genes for enrichment
total_g <- 20501 #total number of genes in the data
path_g <- as.numeric(nrow(features)) # number of genes in a specific pathway
sig_g <- as.numeric(nrow(all_features_results)) #number of total significant prognostic genes in a specific cancer 
sel_g <- as.numeric(nrow(sel_features_results)[1]) ##number of total significant genes found common in pathway

sig_g1 <- sig_g - sel_g #no .pf significant genes which are not in specific pathway
path_g1 <- path_g - sel_g #pathway genes which are not in significant
rest_g <-  path_g1 + sel_g + sig_g1 #significant and genes of pathway
total2_g <- total_g - rest_g  ## genes which are not significant and not in pathway
perc_path <- (sel_g/path_g) * 100
perc_path


#prepare contigenecy table
contingencyTable <- data.frame(
  Sig_Yes=c(sel_g, sig_g1),
  Sig_No=c(path_g1, total2_g))

#row names
row.names(contingencyTable) <- c("PathwayYes", "PathwayNO")
#print contigenecy table
contingencyTable

#apply fisher exact test
Fisher_res <- fisher.test(contingencyTable, alternative = "greater") #one sided
pval_f <- Fisher_res$p.value
odd_ratio <- Fisher_res$estimate

#create dataframe considering important things in fisher test
path_Fisher_test <- cbind(path_g, sel_g, perc_path, pval_f, odd_ratio[1])
#colnames
#colnames(path_Fisher_test) <- c("Total_path_genes" , "Path_sig_Yes", "Percentage", "p.value", "odd_ratio")
#row names
row.names(path_Fisher_test) <- NULL
path_Fisher_test

######################################################
  write.table(cbind(paste0(subfolder[j]),coeff, HR1, CI, pval,path_Fisher_test),
              file=paste0(path,folder[c], "/", "Pos_Pathways_PI_Survival_results.txt"),row.names=F,col.names=F,sep = '\t',append = T);#output file



### Glmnet model  for Lamda
#cross validation model
cvfit1 <- cv.glmnet(as.matrix(train_feature_mat1[3:ncol(train_feature_mat1)]), surv_object_tr, family = "cox", type.measure = "C", nfolds = 3)
lambda_min <- cvfit1$lambda.min

### COX model plot
jpeg(file="Cox_Regression_lamda_plot.jpeg", units="in", width=10, height=10, res=300)
plot(cvfit1)
dev.off()

#Create ROC plot for 1-,3- and 5-years survival time prediction
  tr_roc1 <- survivalROC(Stime        = LUSC_C_tr$OS.time,
                         status       = LUSC_C_tr$OS,
                         marker       = tr_PI$PI,
                         predict.time = 365,
                         method       = "NNE", lambda = lambda_min , span = NULL, window ="symmetric")


  tr_roc3 <- survivalROC(Stime        = LUSC_C_tr$OS.time,
                        status       = LUSC_C_tr$OS,
                         marker       = tr_PI$PI,
                         predict.time = 1095,
                         method       = "NNE", lambda = lambda_min, span = NULL, window ="symmetric")

  tr_roc5 <- survivalROC(Stime        = LUSC_C_tr$OS.time,
                        status       = LUSC_C_tr$OS,
                         marker       = tr_PI$PI,
                         predict.time = 1825,
                         method       = "NNE", lambda = lambda_min, span = NULL, window ="symmetric")

  jpeg(file="ROC.jpeg", units="in", width=10, height=10, res=300)
  plot(tr_roc1$FP, tr_roc1$TP, type="l", xlim=c(0,1), ylim=c(0,1),
       xlab="FP", col="red",
       ylab="TP",main= "AUC Curve for Survival Prediction")
  lines(tr_roc3$FP, tr_roc3$TP, type="l", lty=2, col="blue")
  lines(tr_roc5$FP, tr_roc5$TP, type="l", lty=2, col="green")
  legend(0.5,0.5, legend=c( paste( "1 Year AUC = ",round(tr_roc1$AUC,2)),  paste("3 Years AUC = ",round(tr_roc3$AUC,2)),  paste("5 Years AUC = ",round(tr_roc5$AUC,2))), col =c ("red","blue", "green"), lty=c(1,2), bty="n")

  abline(0,1)
  dev.off()




}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})  
  
}


}
