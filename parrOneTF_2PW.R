#parelle implementation
library(spls)
library(plyr)
library(foreach)
library(doParallel)
library(infotheo)
#library(mpmi)
source("./load_background_data.R")
source("./analyse_functions.R")

y1_y2_pair_alpha = 0.005

################################################################################
#Load Data and convert to data frames as samples in rows and genes in columns
################################################################################
# location of file 1: Gene Expression dataset
file1 <- "./pathways/pws_in_stem/exp.data/At_stem_rma_july2012_AGI129_final.txt"
#file1 <- "./pathways/pws_in_leafs/exp.data/Dataset_drought_stress_seq17_2016st_136samples.csv"
# file 2: Pathway gene IDs
file2 <- "pathways/all_pw/AT_all_pwg.csv"

#file2 <- "./pathways/pws_in_leafs/anthocynin/anthocy_pathway_genes_136.csv"

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
#write.csv(pw.exp,"pw.lignin.csv")
#extract expression for the seleted pathway
pw.exp<-merge(AT_pwlist,exp.data,by = "Gene")
pw.exp<-data.frame(merge(pw.exp,anno,by="Gene"))
pw.exp <- unique(pw.exp)
pw.exp <- within(pw.exp, rm(Gene))
genes <- pw.exp$sym
pw.exp <- as.data.frame(t(pw.exp[,-(sample_size+1)]))#change this
colnames(pw.exp) <- genes
colnames(pw.exp)

################################################################################
#Dimension Reduction Step (Optional)
################################################################################

#obtain spls coefficients
#beta <- spls.coeff(tf.exp,pw.exp)

#t-statisic for spls coefficients
#boot.tstat <- coeff.boot(tf.exp,pw.exp)

#Selected_TF<-data.frame(sort(table(unlist(beta[1:100,seq(1,dim(beta)[2],2)], use.names=FALSE)),decreasing = T))
#Selected_TF<-data.frame(Selected_TF[Selected_TF[,2]>2,])
#Selected_TF <- as.character(Selected_TF$Var1)
#Selected_TF <- Selected_TF[1:5]
#commnet this line to used Dimention reduciton
Selected_TF <- colnames(tf.exp)

#############################################
#code to generate the combination grid of 2 TF and 1 PW gene
pw_pair_tf_all<-create_empty_table(0,11)
colnames(pw_pair_tf_all)<-c("X1","X2","X3","X4","X5","X6","X7","X8","X9", "X10","X11")
for (tf in unique(Selected_TF)){

PW_pairs_df<-data.frame(expand.grid.unique(unique(colnames(pw.exp)),unique(colnames(pw.exp))))#helper function expand.grid.unique is used.
temp_df<-cbind(data.frame(PW_pairs_df,rep(tf,dim(PW_pairs_df)[1])))
colnames(temp_df)<-c("pw1","pw2","tf")

cl<-makeCluster(16)
registerDoParallel(cl)
ls<-foreach(i=1:nrow(temp_df),.combine =rbind,.packages=c("infotheo")) %dopar% {
  vec<-temp_df[i,]
  pw_2tf_1_calc(vec)
  
}
stopCluster(cl)

ls <- data.frame(ls)
rownames(ls)<-NULL
pw_pair_tf_all <-rbind(pw_pair_tf_all,ls)
}

write.csv(pw_pair_tf_all, "PW2TF1_interaction.csv")
