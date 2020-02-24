library(tidyverse)
mytheme <-  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  # theme(panel.border = element_blank(), #Í¼Æ¬Ã»ÓÐ±ß¿ò
  #       axis.line.y = element_line(colour="black",size=0.5, lineend="square"),
  #       axis.line.x = element_line(colour="black",size=0.5, lineend="square"))+
  theme(axis.text.x = element_text(colour="black", size=8))+
  theme(axis.text.y = element_text(colour="black", size=8))+
  theme(axis.title=element_text(size=8))+ theme(legend.position="none")+
  theme(plot.title=element_text(size=8, lineheight=.9))+
  theme(axis.ticks=element_line(colour="black",size=0.5),axis.ticks.length=unit(.2,"lines"))
setwd("E:/ChenYing/postdoctor/²ÄÁÏ/")
# input race-biased gene identified in GTEx
d <- read.table("Data_LUADRaceBiasedGene.tsv",sep="\t",header=T)
head(d)
names(d) <- c("sample","samples","type","age","race","gender",names(d)[7:length(names(d))])
# count sample size of each race
library(tidyverse)
d %>% group_by(race,type) %>% summarise(num=length(samples))
# convert wide to length
library(reshape2)
d_con <- melt(d, id.vars=c("sample","samples","type","age","race","gender"),
              variable.name="name", value.name="fpkm_log")
head(d_con)
# calaulate p value
gene <- names(d)[7:length(names(d))]
p <- data.frame(matrix(c(1:(length(gene)*8)),ncol=8))
head(p)
names(p) <- c("name","pasvsaf","FCasvsaf",
              "pasvseu","FCasvseu",
              "pafvseu","FCafvseu","pnorvstum")
for (i in 1:length(gene)){
  asian <- d_con$fpkm_log[d_con$race=="asian" & d_con$name==gene[i]]
  african <- d_con$fpkm_log[d_con$race=="black or african american" & d_con$name==gene[i]]
  european <- d_con$fpkm_log[d_con$race=="white" & d_con$name==gene[i]]
  p$pasvsaf[i] <- wilcox.test(asian,african)$p.value
  p$FCasvsaf[i] <- mean(african,na.rm=T)-mean(asian,na.rm=T)
  p$pasvseu[i] <- wilcox.test(asian,european)$p.value
  p$FCasvseu[i] <- mean(european,na.rm=T)-mean(asian,na.rm=T)
  p$pafvseu[i] <- wilcox.test(african,european)$p.value
  p$FCafvseu[i] <- mean(european,na.rm=T)-mean(african,na.rm=T)
  p$pnorvstum[i] <- wilcox.test(d_con$fpkm_log[d_con$type=="Solid Tissue Normal" & d_con$name==gene[i]],
                                d_con$fpkm_log[d_con$type=="Primary Tumor" & d_con$name==gene[i]])$p.value
}
p$name <- gene
head(p)
d_con[d_con$name==gene[i],-c(1,2,4,5,6)]
head(d_con)
## plot
rank <- 0
library("grid")
pdf("20200220_luad_race_plot2.pdf",useDingbats = F)
for(i in 1:length(gene))
{
  if(rank%%16==0){
    grid.newpage()
    rank <- 0
  }
  rank <- rank +1
  print(
    ggplot()+geom_jitter(data= subset(d_con,name==gene[i])[,c("race","fpkm_log")],
      aes(race,fpkm_log,colour=race),size=0.2,height = 0)+mytheme+
      ggtitle(gene[i])+ scale_x_discrete(limits=
                                           c("black or african american",
                                             "white"),label=c("African","European"))+
      annotate("text",  x=-Inf, y=Inf, hjust=-.2, vjust=2,
               label=signif(p$pafvseu[p$name==gene[i]],digits=2), vjust=3, size=3)+
      ylab("log2(FPKM-UQ+1)")+#scale_colour_manual(values=c("#FF7F00", "#104E8B"))+
      xlab("Race"),
    vp=viewport(.19,.21,x=(rank %% 4)*0.25+0.12,y=1-(rank %/% 4.1+0.5)*0.25)
  )
}
dev.off()
