setwd("~/projects/ACE2/FPKM-UQ_processTCGA2/ACE2_expr_cor/")
# input ACE2 median expression in 19 cancer type (downloaded from UCSC Xena)
expr <- read.table("Data_MedianExpressionInDifferentCancerType.txt",sep="\t",header=T)
head(expr)
expr1 <- na.omit(expr)
# include all cancer types 
p1 <- ggplot(data=expr,aes(y=Tumor,x=Normal))+geom_point()+mytheme_violin+ geom_abline(intercept=0, slope=1)+
  xlab("Expression level in normal tissue (log2(FPKM+1))")+ylab("Expression level in tumor (log2(FPKM+1))")+
  xlim(10,20)+ylim(10,20)+
  annotate("text", x=-Inf, y=Inf, 
           label=paste("R = ",signif(cor.test(expr$Tumor,expr$Normal)$estimate,digit=2),sep=""), hjust=-.2, vjust=2)+
  annotate("text", x=-Inf, y=Inf, 
           label=paste("P = ",signif(cor.test(expr$Tumor,expr$Normal)$p.value,digit=2),sep=""), hjust=-.2, vjust=4)
cor.test(expr$Tumor,expr$Normal)  
# only include cancer types of which the sample sizes of both primary tumors and solid normal tissues are greater than three
p2 <- ggplot()+geom_point(data=subset(expr,Normal_num > 3),aes(y=Tumor,x=Normal))+mytheme_violin+ geom_abline(intercept=0, slope=1)+
  #xlab("Expression level in normal tissue (log2(FPKM+1))")+ylab("Expression level in tumor (log2(FPKM+1))")+
  xlim(9,20)+ylim(9,20)+
  geom_point(data=expr[19,],aes(y=Tumor,x=Normal),colour="red")+
  geom_point(data=expr[18,],aes(y=Tumor,x=Normal),colour="red")+
  annotate("text", x=-Inf, y=Inf, 
           label=paste("R = ",signif(cor.test(subset(expr,Normal_num > 3)$Tumor,subset(expr,Normal_num > 3)$Normal)$estimate,digit=2),sep=""), hjust=-.2, vjust=2)+
  annotate("text", x=-Inf, y=Inf, 
           label=paste("P = ",signif(cor.test(subset(expr,Normal_num > 3)$Tumor,subset(expr,Normal_num > 3)$Normal)$p.value,digit=2),sep=""), hjust=-.2, vjust=4)#+
  #geom_point(data=subset(expr,Type %in% "LUSC"|Type %in% "LUAD"),aes(y=Tumor,x=Normal),colour="red")
p2
length(subset(expr,Normal_num > 3)$Type)
library("grid")
pdf("ACE2_expr_cor_across_cancer2.pdf",width=20,height=20,useDingbats = F)
print(p1,vp=viewport(.09,.09,x=.2,y=.9))
print(p2,vp=viewport(.09,.09,x=.34,y=.9))
dev.off()
