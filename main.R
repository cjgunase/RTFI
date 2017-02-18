library(spls)
library(plyr)
library(foreach)
library(doParallel)
library(infotheo)
require(fitdistrplus)
#library(mpmi)
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
               "MYB86","VND4","HB53","MYB52","XND1","MYB92")
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

###Read directly##
pw.exp <- read.csv("./pathways/weng_ping_arabi/pathway_data_data1.csv")
tf.exp <- read.csv("./pathways/weng_ping_arabi/tf_data_data2.csv")

posTFs <- as.character(colnames(tf.exp)[1:35])
################################################################################
#Dimension Reduction Step (Optional)
################################################################################

#obtain spls coefficients
beta <- spls.coeff(tf.exp,pw.exp)

#t-statisic for spls coefficients
#boot.tstat <- coeff.boot(tf.exp,pw.exp)

Selected_TF<-data.frame(sort(table(unlist(beta[1:100,seq(1,dim(beta)[2],2)], use.names=FALSE)),decreasing = T))
Selected_TF<-data.frame(Selected_TF[Selected_TF[,2]>2,])
Selected_TF <- Selected_TF$Var1
#commnet this line to used Dimention reduciton
#Selected_TF <- colnames(tf.exp)

################################################################################
#CMI calculations
################################################################################
n=128

vec<-c()
for(i in 1:1000){
  vec<-c(vec,mi3(runif(n),runif(n),runif(n)))
}
hist(vec,main="3Gene interaction for randomised data",xlab="S7/S1+S2+S3")
fit<-fitdist(vec,"norm",method = c("mle"))
plot(fit)

count_i <-0
  for(i in 1:(dim(pw.exp)[2]-1)){
    for(j in (i+1):(dim(pw.exp)[2])){
      for(k in Selected_TF){
      #for(k in j:10){ #comment this for server
        count_i <-count_i + 1
        print(count_i)
        y1 <- pw.exp[,i]
        y2 <- pw.exp[,j]
        x <- tf.exp[[k]]
        #Mutual Information
        y1 <- discretize(y1)
        y2 <- discretize(y2)
        x <- discretize(x)
        
        enty1 <- entropy(y1)
        if(enty1==0) next
        enty2 <- entropy(y2)
        if(enty2==0) next
        entx <- entropy(x)
        if(entx==0) next
        
        enty1x <- condentropy(y1,x, method="emp")
        if(enty1x==0) next
        
        enty2x <- condentropy(y2,x, method="emp")
        if(enty2x==0) next
        
        enty1_y2x <- condentropy(y1,data.frame(x,y2), method="emp")
        if(enty1_y2x==0) next
        
        enty2_y1x <- condentropy(y2,data.frame(x,y1), method="emp")
        if(enty2_y1x==0) next
        
        entx_y1y2 <- condentropy(x,data.frame(y2,y1), method="emp")
        if(entx_y1y2==0) next
       
        #calculates MI y1;y2
        Iy1y2 <- condinformation(y1, y2, method="emp")
        if(Iy1y2==0) next

        #calculates pvalue for MI y1;y2
                
        #calculates MI y1;y2|x
        Iy1y2_x <- condinformation(y1, y2,x, method="emp")
        if(Iy1y2_x==0) next
        
        Iy1x_y2<- condinformation(y1,x,y2, method="emp")
        if(Iy1x_y2==0) next
        
        Iy2x_y1 <- condinformation(y2,x,y1, method="emp")
        if(Iy2x_y1==0) next
        
        Iy1y2x <- Iy1y2 - Iy1y2_x 
        if(Iy1y2x <= 0) next
        I3<-Iy1y2x/(enty1_y2x+enty2_y1x+entx_y1y2)
        I3 <- 0.1176310
        z<-(I3-mean(vec))/sd(vec)
        pnorm(z,lower.tail=FALSE)
        I3_pval <- pnorm(z,lower.tail=FALSE)
        
        print(paste(c(colnames(pw.exp)[i],
                    colnames(pw.exp)[j],
                    k,
                    enty1,
                    enty2,
                    entx,
                    enty1_y2x,
                    enty2_y1x,
                    entx_y1y2,
                    Iy1x_y2,
                    Iy2x_y1,
                    Iy1y2_x,
                    Iy1y2x,
                    I3,
                    I3_pval),
                    sep = ","
        ))
        
        cat(
          	colnames(pw.exp)[i],
          	colnames(pw.exp)[j],
          	k,
	        enty1,
	        enty2,
	        entx,
	        enty1_y2x,
	        enty2_y1x,
	        entx_y1y2,
	        Iy1x_y2,
	        Iy2x_y1,
	        Iy1y2_x,
	        Iy1y2x,
	        I3,
	        I3_pval,
          "\n",
          file="./all_combination.txt",
          sep="\t",
          append=TRUE)
      }
    }
  }

####### replace the file from the server ###############
results <- read.csv("./all_combination.txt",sep="\t",header = F)
results$V16 <- NULL
colnames(results)<-c("Y1","Y2","X","H_Y1","H_Y2","H_X",
                     "S1","S2","S3","S4","S5","S6","S7","S7_div_123","S7_div_123_pval")
########################################################
results$S7_div_sum_1_2_3 <- results$S7/(results$S1+results$S2+results$S3)
hist(results$S7_div_sum_1_2_3)
fit<-fitdist(results$S7_div_sum_1_2_3,"exp",method = c("mle"))
plot(fit)
summary(fit)[1]
crit_S7_div_sum_1_2_3<-qexp(0.95,summary(fit)[1]$estimate)

results<-results[with(results, order(S7_div_123_pval)), ]
hist(results$S7_div_123_pval)
corr.pval<-p.adjust(results$S7_div_123_pval, method = "bonferroni", n = length(results$S7_div_123_pval))
results$corrected.pval<-corr.pval

selected_results <-results[results$corrected.pval < 0.05,]
tf_freq<-data.frame(sort(table(selected_results$X),decreasing = T))
tf_freq<-tf_freq[tf_freq$Freq>0,]


pair_list<-rbind(data.frame(tf=selected_results$X,pw=selected_results$Y1,str=selected_results$S7_div_sum_1_2_3),
                 data.frame(tf=selected_results$X,pw=selected_results$Y2,str=selected_results$S7_div_sum_1_2_3))


#aggregate pair_list df to get average str for each pair
pair_list <- aggregate(pair_list[,3], list(pair_list$tf,pair_list$pw), mean)
colnames(pair_list) <-c("tf","pw","str")


barplot(sort(table(as.character(pair_list$tf)),decreasing = T)[1:20],ylab="pairs frequency",xlab="",
        las=2,ylim = c(0,25),main="Rank TFs based on Pairs Frequency")
TF_rank_MI<-data.frame(sort(table(as.character(pair_list$tf)),decreasing = T))

write.csv(pair_list,"STEM_128.csv")

########### ROC evaluation ########

# 3GInteraction method.
cutoffs<-seq(1,dim(tf_freq)[1],1)
topk<-c()
noPos<-c()
tpr_list<-c()
fpr_list<-c()
for(cutoff in cutoffs){
  k = cutoff
  P<-posTFs
  N<-setdiff(tf_freq$Var1,P)
  
  topK<-as.character(tf_freq[,"Var1"][1:k])
  bottomK<-as.character(tf_freq[,"Var1"][(k+1):dim(tf_freq)[1]])
  
  TP<-intersect(topK,P)
  TPR<-length(TP)/length(P)
  TN<-intersect(N,bottomK)
  SPC<-length(TN)/length(N)
  tpr_list<-c(tpr_list,TPR)
  fpr_list<-c(fpr_list,1-SPC)
}

plot(fpr_list,tpr_list,type="l",col="red",
     xlab="False Positive Rate(FPR)",ylab="True Positive Rate(TPR)",main="ROC Lignin Biosynthesis Pathway",lwd=1.5)

#Embryonic Stem Cells Atlas of Pluripotency Evidence
#Bottom up GGM Partial Correlation
library(stats)
library(MASS)
sig=0.05         ## significance level to test the correlation and partial correlation , can be changed
sig1=0.05        ## significance level to test the difference of correlations, can be changed
sig2=0.05        ##significance level for cutt off intf_tf--------- 
sig_diff=0.05    ## for dimention reduction, triplet with correlation difference less than diff_cutoff are discarded
source("./pcor_func.R")


pwg<-colnames(pw.exp)
tfg<-as.character(Selected_TF)
pw<-as.matrix(t(pw.exp))
tf<-tf.exp[tfg]
tf<-as.matrix(t(tf))
results<-approach1(pw,tf,pwg,tfg)
results<- results[results$regulator.no_regulator=="regulator",]
results$tf
TF_rank_pCor<-data.frame(sort(table(as.character(results$tf)),decreasing = T))


TF_rank_pCor<-TF_rank_pCor[1:dim(TF_rank_MI)[1],]

cutoffs<-seq(1,dim(TF_rank_pCor)[1],1)
topk<-c()
noPos<-c()
tpr_list.2<-c()
fpr_list.2<-c()
for(cutoff in cutoffs){
  k = cutoff
  P<-posTFs
  N<-setdiff(TF_rank_pCor$Var1,P)
  
  topK<-as.character(TF_rank_pCor[,"Var1"][1:k])
  bottomK<-as.character(TF_rank_pCor[,"Var1"][(k+1):dim(TF_rank_pCor)[1]])
  
  TP<-intersect(topK,P)
  TPR<-length(TP)/length(P)
  TN<-intersect(N,bottomK)
  SPC<-length(TN)/length(N)
  tpr_list.2<-c(tpr_list.2,TPR)
  fpr_list.2<-c(fpr_list.2,1-SPC)
}
lines(fpr_list.2,tpr_list.2,col="green",type = "l",lwd=1.5)
legend(x="bottomright",c("3GMI","BottomUpGGM"),lty = c(1,1),lwd=c(2.5,2.5),col=c("red","green"))








