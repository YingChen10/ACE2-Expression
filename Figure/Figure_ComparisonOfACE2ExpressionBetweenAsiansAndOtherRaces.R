library("tidyverse")
mytheme <-  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(colour="black", size=8))+
  theme(axis.text.y = element_text(colour="black", size=8))+
  theme(axis.title=element_text(size=8))+
  #theme(plot.title=element_text(size=1.5, lineheight=.9))+
  theme(axis.ticks=element_line(colour="black",size=0.25),axis.ticks.length=unit(.2,"lines"))+
  theme(legend.title=element_text(size=8))+
  theme(legend.text=element_text(size=8))
setwd("~/projects/ACE2/cancer_type_7/")
# input ACE2 expression level in 7 cancer type (log2(FPKM_UQ+1)); downloaded from UCSC Xena
type <- read.table("Data_ACE2_expression_among_seven_cancer_type.txt",sep="\t",header=T)
head(type)
names(type) <- c("type","sample","samples","tissue","age","gender","race","expr")
type$expr <- as.numeric(as.character(type$expr))
head(type)
str(type)
## other races marked with NA
type$race[type$race == "not reported" |type$race == "american indian or alaska native" |
            type$race == "native hawaiian or other pacific islander"] <- NA
## ANOVA 
summary(a <- lm(type$expr ~ type$type + type$age + type$gender + type$race))
anova(a)
# subset type
type1 <- type[,c("type","sample","race","expr")]
head(type1)
type2 <- na.omit(type1)
unique(type2$race)
## combine Aferican and European and mark with "Others"
type2$race1 <- "Others"
type2$race1[type2$race=="asian"] <- "Asian"

##compare ACE2 expression between Asians and others and calculate the p value 
## calculate the sample size of each race for each cancer type
subset(type2,race1=="Asian")[,sort("expr")]
a <- as.data.frame(type2 %>% group_by(type,race1) %>% summarise(median=median(expr)))
p <- vector(length=length(unique(type2$type)))
asian <- vector(length=length(unique(type2$type)))
other <- vector(length=length(unique(type2$type)))
t <- unique(type2$type)
for(i in 1:length(unique(type2$type))){
  p[i] <- wilcox.test(subset(type2,type==t[i] & race1=="Asian")$expr,subset(type2,type==t[i] & race1=="Others")$expr)$p.value
  asian[i] <- length(subset(type2,type==t[i] & race1=="Asian")$expr)
  other[i] <- length(subset(type2,type==t[i] & race1=="Others")$expr)
}
info <- data.frame(type=unique(type2$type),p=p,asian=asian,other=other)
info

# plot 
p1 <- ggplot(data=type2,aes(x=reorder(type,expr),y=expr,colour=race1))+
  geom_boxplot(position="dodge",outlier.size = 0.2,outlier.colour = "white")+
  mytheme+ylim(0,25)+
  scale_colour_manual(values=c("#EE7600", "#104E8B"))+
  #scale_x_discrete(limits=info$type)
  annotate("text", x=1, y=22, label=paste("P = ",signif(info$p[1],digit=1),sep="")) +
  annotate("text", x=2, y=22, label=paste("P = ",signif(info$p[2],digit=1),sep=""))+
  annotate("text", x=3, y=22, label=paste("P = ",signif(info$p[3],digit=1),sep=""))+
  annotate("text", x=4, y=22, label=paste("P = ",signif(info$p[4],digit=1),sep=""))+
  annotate("text", x=5, y=22, label=paste("P = ",signif(info$p[5],digit=1),sep=""))+
  annotate("text", x=6, y=22, label=paste("P = ",signif(info$p[6],digit=1),sep=""))+
  annotate("text", x=7, y=22, label=paste("P = ",signif(info$p[7],digit=1),sep=""))+
  annotate("text", x=8, y=22, label=paste("P = ",signif(info$p[8],digit=1),sep=""))+
  annotate("text", x=9, y=22, label=paste("P = ",signif(info$p[9],digit=1),sep=""))+
  annotate("text", x=1, y=5, label=info$asian[1])+
  annotate("text", x=2, y=5, label=info$asian[2])+
  annotate("text", x=3, y=5, label=info$asian[3])+
  annotate("text", x=4, y=5, label=info$asian[4])+
  annotate("text", x=5, y=5, label=info$asian[5])+
  annotate("text", x=6, y=5, label=info$asian[6])+
  annotate("text", x=7, y=5, label=info$asian[7])+
  annotate("text", x=8, y=5, label=info$asian[8])+
  annotate("text", x=9, y=5, label=info$asian[9])+
  annotate("text", x=1, y=4, label=info$other[1])+
  annotate("text", x=2, y=4, label=info$other[2])+
  annotate("text", x=3, y=4, label=info$other[3])+
  annotate("text", x=4, y=4, label=info$other[4])+
  annotate("text", x=5, y=4, label=info$other[5])+
  annotate("text", x=6, y=4, label=info$other[6])+
  annotate("text", x=7, y=4, label=info$other[7])+
  annotate("text", x=8, y=4, label=info$other[8])+
  annotate("text", x=9, y=4, label=info$other[9])
p1
library("grid")
pdf("20200207FigS5.pdf",width=20,height=20,useDingbats = F)
print(p1,vp=viewport(.2,.1,x=.2,y=.9))
dev.off()

  




