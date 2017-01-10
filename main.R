library(spls)
library(plyr)
library(foreach)
library(doParallel)
library(infotheo)
library(mpmi)
source("./load_background_data.R")
source("./analyse_functions.R")

alpha = 0.05

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
#beta <- spls.coeff(tf.exp,pw.exp)

#t-statisic for spls coefficients
#boot.tstat <- coeff.boot(tf.exp,pw.exp)

#Selected_TF<-data.frame(sort(table(unlist(beta[1:150,seq(1,dim(beta)[2],2)], use.names=FALSE)),decreasing = T))
#Selected_TF<-data.frame(test[test[,2]>2,])
#Selected_TF <- rownames(Selected_TF)
#commnet this line to used Dimention reduciton
Selected_TF <- colnames(tf.exp)

################################################################################
#CMI calculations
################################################################################

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
        Iy1y2_pval <- 2*pnorm(-abs(cmi.pw(y1,y2)$zvalue))
        
        #calculates MI y1;y2|x
        Iy1y2_x <- condinformation(y1, y2,x, method="emp")
        if(Iy1y2_x==0) next
        
        Iy1x_y2<- condinformation(y1,x,y2, method="emp")
        if(Iy1x_y2==0) next
        
        Iy2x_y1 <- condinformation(y2,x,y1, method="emp")
        if(Iy2x_y1==0) next
        
        Iy1y2x <- Iy1y2 - Iy1y2_x 
        if(Iy1y2x <= 0) next
        
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
                    Iy1y2,
                    Iy1y2_pval,
                    Iy1y2x),
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
	        Iy1y2,
	        Iy1y2_pval,
	        Iy1y2x,
          "\n",
          file="./all_combination_calculations.txt",
          sep="\t",
          append=TRUE)
      }
    }
  }

####### replace the file from the server ###############
results <- read.csv("./SPLS_combination_calculations_stem_lignin_pw.txt",sep="\t",header = F)
########################################################
#correction for p_value for Iy1_y2_pval

corrected.p.value<-p.adjust(results$V14, method = "BH", n = length(results$V14))
results$V16 <- corrected.p.value
colnames(results)<-c("Y1","Y2","X","H_Y1","H_Y2","H_X",
                     "S1","S2","S3","S4","S5","S6","S6_S7","S6_S7_pval","S7","S6_S7_adj.pval")

colnames(results)

results$S7_div_sum_1_2_3 <- results$S7/(results$S1+results$S2+results$S3)
hist(results$S7_div_sum_1_2_3)
fit<-fitdist(results$S7_div_sum_1_2_3,"exp",method = c("mle"))
plot(fit)
summary(fit)[1]
crit_S7_div_sum_1_2_3<-qexp(0.95,summary(fit)[1]$estimate)

results$S4_div_sum_1_4_6_7 <- results$S4/(results$S1+results$S4+results$S6+results$S7)
hist(results$S4_div_sum_1_4_6_7)
fit<-fitdist(results$S4_div_sum_1_4_6_7,"lnorm",method = c("mle"))
plot(fit)
summary(fit)[1]
#crit_S4_div_sum_1_4_6_7 <- qnorm(0.95,mean = summary(fit)[1]$estimate[1],sd=summary(fit)[1]$estimate[2],lower.tail = T)
crit_S4_div_sum_1_4_6_7 <- qlnorm(0.95,meanlog = summary(fit)[1]$estimate[1],sdlog=summary(fit)[1]$estimate[2],lower.tail = T)

results$S5_div_sum_2_5_6_7 <- results$S5/(results$S2+results$S5+results$S6+results$S7)
hist(results$S5_div_sum_2_5_6_7)
fit<-fitdist(results$S5_div_sum_2_5_6_7,"lnorm",method = c("mle"))
plot(fit)
summary(fit)[1]
#crit_S5_div_sum_2_5_6_7 <- qnorm(0.95,mean = summary(fit)[1]$estimate[1],sd=summary(fit)[1]$estimate[2],lower.tail = T)
crit_S5_div_sum_2_5_6_7 <- qlnorm(0.95,meanlog = summary(fit)[1]$estimate[1],sdlog=summary(fit)[1]$estimate[2],lower.tail = T)


























