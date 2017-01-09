library(spls)
library(plyr)
library(foreach)
library(doParallel)
library(infotheo)
library(mpmi)
source("./load_background_data.R")
source("./analyse_functions.R")

################################################################################
#Load Data and convert to data frames as samples in rows and genes in columns
################################################################################
# location of file 1: Gene Expression dataset
file1 <- "./pathways/pws_in_stem/exp.data/At_stem_rma_july2012_AGI129_final.txt"
#file1 <- "./pathways/pws_in_leafs/exp.data/Dataset_drought_stress_seq17_2016st_136samples.csv"
# file 2: Pathway gene IDs
file2 <- "./pathways/pws_in_stem/lignin/lignin_pathway_genes.csv"
#file2 <- "./pathways/pws_in_leafs/anthocynin/anthocy_pathway_genes.csv"

#todo:
#extract the directory from file2 and put resutls in the same directory

# file3 : known positive TFs for this pathway
#file3 <- "./pathways/pws_in_leafs/anthocynin/anthocy_pTFs.csv"
#load user data
##############################################################################
#1. User need to load gene expression dataset
exp.data<-read.csv(file1,sep = "\t",stringsAsFactors=FALSE,header = T)
sample_size<-dim(exp.data)[2]-1
#2. User loads Interested Pathway
AT_pwlist <-read.csv(file2)
AT_pwlist<-data.frame(AT_pwlist$pw)
colnames(AT_pwlist)<-"Gene"

#3.if known pTFs(upload in positiveTF symbols)
#pTFlist <-read.csv(file3,header = T)
#posTFs<-data.frame(pTFlist$sym)
#posTFs<-as.vector(posTFs[,1])

#users can specify known TFs which regulates the interested pathway

posTFs<-c("GATA12","LBD15","MYB103","MYB43","MYB63","SND3","MYB46","MYB61",
               "LBD30","MYB85","SND1","NST1","MYB58","NST2","SND2","GRF3","SHP1",
               "MYB86","VND4","HB53","MYB52","XND1")

##############################################################################

#extract expression for the TFs
tf.exp<-merge(AT_TF.list,exp.data,by = "Gene")
tf.exp<-as.data.frame(merge(tf.exp,anno,by="Gene"))
tf.exp<-tf.exp[!duplicated(tf.exp$sym),]
tf.exp <- unique(tf.exp)
tf.exp <- within(tf.exp, rm(Gene))
genes <- tf.exp$sym
tf.exp <- data.frame(t(tf.exp[,-(sample_size+1)]))#change this 129 for stem 137 for leaf
colnames(tf.exp) <- genes
write.csv(pw.exp,"pw.lignin.csv")
#extract expression for the seleted pathway
pw.exp<-merge(AT_pwlist,exp.data,by = "Gene")
pw.exp<-data.frame(merge(pw.exp,anno,by="Gene"))
pw.exp <- unique(pw.exp)
pw.exp <- within(pw.exp, rm(Gene))
genes <- pw.exp$sym
pw.exp <- as.data.frame(t(pw.exp[,-(sample_size+1)]))#change this
colnames(pw.exp) <- genes


################################################################################
#Dimension Reduction Step (Optional)
################################################################################

#obtain spls coefficients
beta <- spls.coeff(tf.exp,pw.exp)

#t-statisic for spls coefficients
#boot.tstat <- coeff.boot(tf.exp,pw.exp)

Selected_TF<-data.frame(sort(table(unlist(beta[1:150,seq(1,dim(beta)[2],2)], use.names=FALSE)),decreasing = T))
Selected_TF<-data.frame(test[test[,2]>2,])



################################################################################
#CMI calculations
################################################################################





























