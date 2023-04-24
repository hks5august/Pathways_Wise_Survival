#Load required libraries
library(caret)
library(ggfortify)
library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(ggfortify)
library(MASS)
args <- commandArgs(TRUE)

#setwd("/Users/kaurh8/Documents/Pathways_Survival_analysis_all_cancer")

#setwd("/Users/kaurh8/Documents/CCGA_datasets/mRNA_693_samples/LGG/LGG/survival")

set.seed(7)

data <- read.table(args[1], header = TRUE, sep = ",", row.names = 1, check.names = FALSE)
#data <- read.table("LIHC.csv", header = TRUE, sep = ",", row.names = 1, check.names = FALSE)
dim(data)

Clin_data <- data[1:14]
head(Clin_data)
dim(Clin_data)
Clin_data1 <- Clin_data %>% mutate(OS_month = round(OS.time/30.417, digit=0))

exp_data <- data[15:ncol(data)]
dim(exp_data)

data1 <- cbind(Clin_data1, exp_data)

# create a file to store results
write.table(cbind("ID","Beta","HR","P-value","GP1","GP2","Hr-Inv-lst","Concordance","Std_Error"),
            file="./Survival_results_for_genes.txt",row.names=F,col.names=F,sep = '\t');

#Here features to compute survival start from 16th column onwards, which is transcriptomics features

for(i in seq(from=16, to=length(data1), by=1))
{
  surv_object2 <- Surv(time = data1$OS_month, event = data1$OS)
  
  #survival analysis: fits cox ph model to find HR for median cut
  fit1 <- survfit(surv_object2 ~ (data1[,i])>(median(data1[1,i])), data=data1);
  summary(fit1);
  fit1.coxph <- coxph(surv_object2 ~ (data1[,i])>(median(data1[1,i])), data = data1)
  # summary(fit1.coxph);
  first <- coef(summary(fit1.coxph))
  
  write.table(cbind(colnames(data1[i]),first[1],first[2],first[5],fit1$n[1],fit1$n[2],1/first[2],fit1.coxph$concordance[6],fit1.coxph$concordance[7]),
               file="./Survival_results_for_genes.txt",row.names=F,col.names=F,sep = '\t',append = T);#output file
 
}


### Significant results 

write.table(cbind("ID","Beta","HR","P-value","GP1","GP2","Hr-Inv-lst","Concordance","Std_Error"),
            file="./Significant_Survival_results_for_genes.txt",row.names=F,col.names=F,sep = '\t');


for(i in seq(from=16, to=length(data), by=1))
{
  surv_object2 <- Surv(time = data1$OS_month, event = data1$OS)
  
  #survival analysis: fits cox ph model to find HR for median cut
  fit1 <- survfit(surv_object2 ~ (data1[,i])>(median(data1[1,i])), data=data1);
  summary(fit1);
  fit1.coxph <- coxph(surv_object2 ~ (data1[,i])>(median(data1[1,i])), data = data1)
  # summary(fit1.coxph);
  first <- coef(summary(fit1.coxph))
  
  #check whether the pvalue is significant (< 0.05) or not
  if((first[5]<=0.05)&&(!is.na(first[5]))&&(!is.na(first[2])))
  {write.table(cbind(colnames(data1[i]),first[1],first[2],first[5],fit1$n[1],fit1$n[2],1/first[2],fit1.coxph$concordance[6],fit1.coxph$concordance[7]),
               file="./Significant_Survival_results_for_genes.txt",row.names=F,col.names=F,sep = '\t',append = T);#output file
  }
}



