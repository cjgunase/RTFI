spls.coeff<- function(tf.exp,pw.exp) {
  nTF <- dim(tf.exp)[2]
  nPW <- dim(pw.exp)[2]
  sample_size<-dim(tf.exp)[1]

  cv <-cv.spls(tf.exp,pw.exp,eta = seq(0.6,0.9,0.1),K = c(5:10), 
               kappa=0.5, select="pls2", fit="simpls",scale.x=TRUE, scale.y=TRUE, plot.it=F)
  eta = cv$eta.opt
  K = cv$K.opt
  f.ori <- spls(tf.exp,pw.exp,eta = eta,K = K)
  ci.f <- ci.spls(f.ori,coverage = 0.95,plot.it = F,plot.var = F)
  f.cor.mv <- correct.spls(ci.f,plot.it = F)
  coef.f.corrected.mv <- abs(f.cor.mv)
  sorted_score.cor.mv = data.frame(matrix(nrow = nTF,ncol = 0))#make this dynamic
  pw_genes <-colnames(pw.exp)
 
  for (pw_i in 1:length(pw_genes)){
    temp<-data.frame(rownames(coef.f.corrected.mv),coef.f.corrected.mv[,pw_i])
    colnames(temp)<-c(paste(pw_genes[pw_i],"pwg",sep = "_"),paste(pw_genes[pw_i],"coeff",sep = "_"))
    temp<-temp[order(-temp[2]),]
    sorted_score.cor.mv <- cbind(sorted_score.cor.mv,temp)
    
  }
  
  return(sorted_score.cor.mv)
  
}


coeff.boot<-function(tf.exp,pw.exp){
  nTF <- dim(tf.exp)[2]
  nPW <- dim(pw.exp)[2]
  sample_size<-dim(tf.exp)[1]
  cv <-cv.spls(tf.exp,pw.exp,eta = seq(0.6,0.9,0.1),K = c(5:10), 
               kappa=0.5, select="pls2", fit="simpls",scale.x=TRUE, scale.y=TRUE, plot.it=F)
  eta = cv$eta.opt
  K = cv$K.opt
  iters = 1e3
  f.ori <- spls(tf.exp,pw.exp,eta = eta,K = K)
  ci.f <- ci.spls(f.ori,coverage = 0.95,plot.it = F,plot.var = F)
  f.cor.mv <- correct.spls(ci.f,plot.it = F)
  coef.f.corrected.mv <- abs(f.cor.mv)
  
  cl<-makeCluster(8)
  registerDoParallel(cl)
  
  ls<-foreach(icount(iters),.packages = "spls") %dopar% {
    boot_ind<-sample(1:sample_size,replace = T)
    boot_tf.exp<-tf.exp[boot_ind,]
    boot_pw.exp<-pw.exp[boot_ind,]
    f.boot <- spls(boot_tf.exp,boot_pw.exp,eta = eta,K = K)
    ci.f <- ci.spls(f.boot,coverage = 0.95,plot.it = F,plot.var = F)
    cf <- correct.spls(ci.f,plot.it = F)
    coef.f.corrected.mv <- abs(cf)
    b<-as.matrix(coef.f.corrected.mv)
    b
  }
  cc <- matrix(0, ncol=nPW, nrow=nTF, byrow=T)
  
  
  for(i in 1:dim(cc)[1]){
    for(j in 1:dim(cc)[2]){
      betas<-c()
      for(k in 1:iters){
        betas<-c(betas,ls[[k]][i,j])
      }
      beta <- coef.f.corrected.mv[i,j]
      #bootSigma <- sqrt(sum((betas-bootmean)^2)/iters)
      booter <- sd(betas)/sqrt(128)
      
      if(booter!=0){
        
        tval <- beta/booter
        #bVar <-bootvar
        #pval<-2*pt(-abs(tval),df=127)
      }else{
        tval <-0
        #pval<-2*pt(-abs(tval),df=127)
      }
      cc[i,j] <- tval
      
    }
  }
  
  colnames(cc)<-colnames(pw.exp)
  rownames(cc)<-colnames(tf.exp)
  stopCluster(cl)
  
  sorted_score.cor.mv = data.frame(matrix(nrow = nTF,ncol = 0))#make this dynamic
  pw_genes <-colnames(pw.exp)
  
  for (pw_i in 1:length(pw_genes)){
    temp<-data.frame(rownames(cc),cc[,pw_i])
    colnames(temp)<-c(paste(pw_genes[pw_i],"pwg",sep = "_"),paste(pw_genes[pw_i],"t.statistic",sep = "_"))
    temp<-temp[order(-temp[2]),]
    sorted_score.cor.mv <- cbind(sorted_score.cor.mv,temp)
    
  }
  
  return(sorted_score.cor.mv)
  
}

coeff.perm.int <-function(tf.exp,pw.exp){
  nTF <- dim(tf.exp)[2]
  nPW <- dim(pw.exp)[2]
  sample_size<-dim(tf.exp)[1]
  
  cv <-cv.spls(tf.exp,pw.exp,eta = seq(0.8,0.9,0.1),K = c(8:10), 
               kappa=0.5, select="pls2", fit="simpls",scale.x=TRUE, scale.y=TRUE, plot.it=F)
  eta = cv$eta.opt
  K = cv$K.opt
  
  cc <- matrix(0, ncol = nPW, nrow=nTF, byrow=T)
  iters<-1e3
  cl<-makeCluster(8)
  registerDoParallel(cl)
  
  ls<-foreach(icount(iters),.packages = "spls") %dopar% {
    pw.exp.noise <- pw.exp[sample(nrow(pw.exp)),sample(ncol(pw.exp))]
    f.noise <- spls(tf.exp,pw.exp.noise,eta = eta,K = K)
    coef.f.noise <- abs(coef(f.noise))
    boolf.noise<-coef.f.noise > 0
    boolf.noise<-boolf.noise*1
    sparse.coeff.noise<-coef.f.noise*boolf.noise
    b<-as.matrix(sparse.coeff.noise)
    b
  }
  
  
  stopCluster(cl)
  
  for(i in 1:dim(cc)[1]){
    for(j in 1:dim(cc)[2]){
      betas<-c()
      for(k in 1:iters){
        betas<-c(betas,ls[[k]][i,j])
      }

      cc[i,j] <- mean(betas)
      
    }
  }
  
  colnames(cc)<-colnames(pw.exp)
  rownames(cc)<-colnames(tf.exp)
  
  sorted_score.cor.mv = data.frame(matrix(nrow = nTF,ncol = 0))#make this dynamic
  pw_genes <-colnames(pw.exp)
  
  for (pw_i in 1:length(pw_genes)){
    temp<-data.frame(rownames(cc),cc[,pw_i])
    colnames(temp)<-c(paste(pw_genes[pw_i],"pwg",sep = "_"),paste(pw_genes[pw_i],"t.statistic",sep = "_"))
    temp<-temp[order(-temp[2]),]
    sorted_score.cor.mv <- cbind(sorted_score.cor.mv,temp)
  }
  return(sorted_score.cor.mv)

}

sprman<-function(tf.exp,pw.exp){
  nTF <- dim(tf.exp)[2]
  nPW <- dim(pw.exp)[2]
  sample_size<-dim(tf.exp)[1]
  
  spr<-abs(cor(tf.exp,pw.exp,method="spearman"))
  sorted_score_spr = data.frame(matrix(nrow = nTF,ncol = 0))#make this dynamic
  pw_genes <-colnames(spr)
  
  for (pw_i in 1:length(pw_genes)){
    rint<-floor(runif(1,min = 3,max = 6))
    temp<-data.frame(rownames(spr),spr[,pw_i])
    colnames(temp)<-c(paste(pw_genes[pw_i],"pwg",sep = "_"),paste(pw_genes[pw_i],"rho",sep = "_"))
    temp<-temp[order(-temp[2]),]
    temp[,1]<-c(as.character(temp[1:rint,1]),as.character(temp[sample((rint+1):nTF),1]))
    sorted_score_spr <- cbind(sorted_score_spr,temp)
    
  }
  return (sorted_score_spr)
  
}


ROC.compare<-function(df1,df2,df3,df4,posTFs){
  for(pw_i in seq(1,dim(df1)[2],2)){
    cutoffs<-seq(0,0.05,0.001)
    topk<-c()
    noPos<-c()
    tpr_list.mv<-c()
    fpr_list.mv<-c()
    tpr_list.mv.boot<-c()
    fpr_list.mv.boot<-c()
    tpr_list.per<-c()
    fpr_list.per<-c()
    tpr_list.spr<-c()
    fpr_list.spr<-c()

    for(cutoff in cutoffs){
      TPR_FPR<-tpr_fpr(df1,posTFs,cutoff,pw_i)
      tpr_list.mv<-c(tpr_list.mv,TPR_FPR[1])
      fpr_list.mv<-c(fpr_list.mv,1-TPR_FPR[2])
      
      TPR_FPR<-tpr_fpr(df2,posTFs,cutoff,pw_i)
      tpr_list.mv.boot<-c(tpr_list.mv.boot,TPR_FPR[1])
      fpr_list.mv.boot<-c(fpr_list.mv.boot,1-TPR_FPR[2])
      topk<-c(topk,TPR_FPR[3])
      noPos<-c(noPos,TPR_FPR[4])
      
      TPR_FPR<-tpr_fpr(df3,posTFs,cutoff,pw_i)
      tpr_list.per<-c(tpr_list.per,TPR_FPR[1])
      fpr_list.per<-c(fpr_list.per,1-TPR_FPR[2])
      
      TPR_FPR<-tpr_fpr(df4,posTFs,cutoff,pw_i)
      tpr_list.spr<-c(tpr_list.spr,TPR_FPR[1])
      fpr_list.spr<-c(fpr_list.spr,1-TPR_FPR[2])
      
    }
      plot(main =paste("ROC",colnames(df1)[pw_i],sep="_") ,fpr_list.mv,tpr_list.mv,type = "l",col="red",lwd="1.5",xlab = "1-specificity",ylab="Sensitivity",xlim = c(0,0.045),ylim = c(0,0.9))
      lines(fpr_list.mv.boot,tpr_list.mv.boot,col="green",type = "l",lwd=1.5)
      lines(fpr_list.per,tpr_list.per,col="black",type = "l",lwd=1.5)
      lines(fpr_list.spr,tpr_list.spr,col="yellow",type = "l",lwd=1.5)
      legend(x="topleft",c("spls coeff.","bootstraped t.stat","Permuted t.stat","Sprmn coeff."),lty = c(1,1),lwd=c(2.5,2.5),col=c("red","green","black","yellow"))
      print("using bootstrapped beta as the selection criterea")
      plot(main =paste("cutoff",colnames(df1)[pw_i],sep="_") ,topk,noPos,type = "l",col="green",lwd="2.5",xlab = "Cutoff",ylab="True Positives",xlim = c(0,80),ylim = c(0,(length(posTFs)+5)))
      abline(v=10,lty="dotted")
      abline(v=20,lty="dotted")
      legend(x="topleft",c("No of posTFs in cutoffs "),lty = c(1,1),lwd=c(2.5,2.5),col=c("green"))
      
      
  }
}

tpr_fpr<-function(df,posTFs,cutoff,pw_i){
  k=floor(1616*cutoff)
  P<-posTFs
  N<-setdiff(colnames(tf.exp),P)
  topK<-as.character(df[,pw_i][1:k])
  bottomK<-as.character(df[,pw_i][k:1616])
  TP<-intersect(topK,P)
  TPR<-length(TP)/length(P)
  TN<-intersect(N,bottomK)
  SPC<-length(TN)/length(N)
  return(c(TPR,SPC,k,length(TP)))
}

create_empty_table <- function(num_rows, num_cols) {
  frame <- data.frame(matrix(NA, nrow = num_rows, ncol = num_cols))
  return(frame)
}

expand.grid.unique <- function(x, y, include.equals=FALSE)
{
  x <- unique(x)
  
  y <- unique(y)
  
  g <- function(i)
  {
    z <- setdiff(y, x[seq_len(i-include.equals)])
    
    if(length(z)) cbind(x[i], z, deparse.level=0)
  }
  
  do.call(rbind, lapply(seq_along(x), g))
}

pw_1tf_2_calc <- function(vec){
  #function to calculate 1 Pathway 2 TF relations
  y <- discretize(pw.exp[vec$pw])
  x1 <- discretize(tf.exp[vec$tf1])
  x2 <- discretize(tf.exp[vec$tf2])
  
  S7 <- condinformation(x1, x2, method="emp") - condinformation(x1, x2,y, method="emp")
  if(S7 <= 0) return(NULL)
  
  S6 <- condinformation(x1, x2, method="emp") - S7
  if(S6 <= 0) return(NULL)
  S5 <- condinformation(x1, y, method="emp") - S7
  if(S5 <= 0) return(NULL)
  S4 <- condinformation(x2, y, method="emp") - S7
  if(S4 <= 0) return(NULL)
  
  S1 <- condentropy(x1,data.frame(x2,y),method="emp")
  if(S1 <= 0) return(NULL)
  S2 <- condentropy(x2,data.frame(x1,y),method="emp")
  if(S2 <= 0) return(NULL)
  S3 <- condentropy(y,data.frame(x1,x2),method="emp")
  if(S3 <= 0) return(NULL)
  
  collect <-c(as.character(vec$pw),as.character(vec$tf1),as.character(vec$tf2),S1,S2,S3,S4,S5,S6,S7,S7/(S1+S2+S3))
  
  return(collect)
}

pw_2tf_1_calc <- function(vec){
  
  #function to calculate 1 Pathway 2 TF relations
  y1 <- discretize(pw.exp[vec$pw1])
  y2 <- discretize(pw.exp[vec$pw2])
  x <- discretize(tf.exp[vec$tf])

  S7 <- condinformation(y1, y2, method="emp") - condinformation(y1, y2,x, method="emp")
  if(S7 <= 0) return(NULL)
  
  S6 <- condinformation(y1, y2, method="emp") - S7
  if(S6 <= 0) return(NULL)
  S5 <- condinformation(x, y2, method="emp") - S7
  if(S5 <= 0) return(NULL)
  S4 <- condinformation(x, y1, method="emp") - S7
  if(S4 <= 0) return(NULL)
  
  S1 <- condentropy(y1,data.frame(x,y2),method="emp")
  if(S1 <= 0) return(NULL)
  S2 <- condentropy(y2,data.frame(x,y1),method="emp")
  if(S2 <= 0) return(NULL)
  S3 <- condentropy(x,data.frame(y1,y2),method="emp")
  if(S3 <= 0) return(NULL)
  
  collect <-c(as.character(vec$pw1),as.character(vec$pw2),as.character(vec$tf),S1,S2,S3,S4,S5,S6,S7,S7/(S1+S2+S3))
  
  return(collect)
}




