library(dplyr)
library(MASS);
args <- commandArgs(TRUE)
set.seed(7)
#setwd("/Users/kaurh8/Documents/CCGA_datasets/mRNA_693_samples/LGG/primary_vs_Recurrent/DGE")


train <- read.table("../full_train_clin_23271.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
test1 <- read.table("../full_test_clin_23271.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
test2 <- read.table("../full_ext_clin_23271.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

dim(train)
dim(test1)
dim(test2)

features <- read.table(args[1], header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
#features <- read.table("Neuroactive_pathway_genes", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
#features <- read.table("Neuroactive_ligand_receptor_interactionpathway_sig_genes.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
#features <- read.table("Pathways_in_cancerpathway_sig_genes.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

#dataset preparation for selected features

dim(features)

sel_ftr <- as.data.frame(row.names(features))

#sel_ftr1 <- rbind(c("Class"), sel_ftr)
sel_ftr1 <- rbind(c("Death_event", "OS_month"), sel_ftr)
names(sel_ftr1) <- c("ID")
head(sel_ftr1)


#Prepare  data with selected features
sel_train<- as.data.frame(train[,colnames(train) %in% c(sel_ftr1$ID), ])
sel_test1<- as.data.frame(test1[,colnames(test1) %in% c(sel_ftr1$ID), ])
sel_test2<- as.data.frame(test2[,colnames(test2) %in% c(sel_ftr1$ID), ])
#sel_data1 <- as.data.frame(LGG_P_vs_R_data2[,colnames(LGG_P_vs_R_data2) %in% c(sel_ftr1$ID), ])

#Check dimensions of dataset
head(sel_train[1:50])

dim(sel_train)
dim(sel_test1)
dim(sel_test2)


sel_train <- na.omit(sel_train)
sel_test1 <- na.omit(sel_test1)
sel_test2 <- na.omit(sel_test2)

dim(sel_train)
dim(sel_test1)
dim(sel_test2)

write.table(cbind("ID"=rownames(sel_train), sel_train),file=args[2],sep="\t",quote=F, row.names=F)
write.table(cbind("ID"=rownames(sel_test1), sel_test1),file=args[3],sep="\t",quote=F, row.names=F)
write.table(cbind("ID"=rownames(sel_test2), sel_test2),file=args[4],sep="\t",quote=F, row.names=F)


#write.table(cbind("ID"=rownames(sel_train), sel_train),file="Neuroactive_train.txt",sep="\t",quote=F, row.names=F)
#write.table(cbind("ID"=rownames(sel_test1), sel_test1),file="Neuroactive_test1.txt",sep="\t",quote=F, row.names=F)
#write.table(cbind("ID"=rownames(sel_test2), sel_test2),file="Neuroactive_test2.txt",sep="\t",quote=F, row.names=F)


#write.table(cbind("ID"=rownames(sel_train), sel_train),file="Cancer_path_train.txt",sep="\t",quote=F, row.names=F)
#write.table(cbind("ID"=rownames(sel_test1), sel_test1),file="Cancer_path_test1.txt",sep="\t",quote=F, row.names=F)
#write.table(cbind("ID"=rownames(sel_test2), sel_test2),file="Cancer_path_test2.txt",sep="\t",quote=F, row.names=F)
