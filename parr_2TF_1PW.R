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
file2 <- "./pathways/pws_in_stem/lignin/lignin_pathway_genes.csv"
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
beta <- spls.coeff(tf.exp,pw.exp)

#t-statisic for spls coefficients
#boot.tstat <- coeff.boot(tf.exp,pw.exp)

Selected_TF<-data.frame(sort(table(unlist(beta[1:50,seq(1,dim(beta)[2],2)], use.names=FALSE)),decreasing = T))
Selected_TF<-data.frame(Selected_TF[Selected_TF[,2]>2,])
Selected_TF <- as.character(Selected_TF$Var1)
#commnet this line to used Dimention reduciton
#Selected_TF <- colnames(tf.exp)

#############################################
#code to generate the combination grid of 2 TF and 1 PW gene
TF_pairs_df<-data.frame(expand.grid.unique(Selected_TF,Selected_TF))#helper function expand.grid.unique is used.
pw_pair_tf_all<-create_empty_table(0,3)
colnames(pw_pair_tf_all)<-c("pw","tf1","tf2")
pw.column_df<-data.frame()
for (pw in unique(colnames(pw.exp))){
  temp_df<-cbind(data.frame(rep(pw,dim(TF_pairs_df)[1]),TF_pairs_df))
  colnames(temp_df)<-c("pw","tf1","tf2")
  pw_pair_tf_all<-rbind(pw_pair_tf_all,temp_df)
}

############################################
#pw_pair_tf_all<-pw_pair_tf_all[sample(1:1310478,replace = F)[1:1000],]

cl<-makeCluster(8)
registerDoParallel(cl)
ls<-foreach(i=1:nrow(pw_pair_tf_all),.combine =rbind,.packages=c("infotheo")) %dopar% {
  vec<-pw_pair_tf_all[i,]
  pw_1tf_2_calc(vec)
  
}
stopCluster(cl)

write.csv(ls, "onePW_2TF_interaction.csv")



  



