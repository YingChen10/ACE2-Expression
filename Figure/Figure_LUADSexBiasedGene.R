setwd("E:/ChenYing/postdoctor/²ÄÁÏ/")
# input sex-biased gene identified in GTEx
d <- read.table("Data_LUADSexBiasedGene.tsv",sep="\t",header=T)
head(d)
names(d) <- c("sample","samples","type","age","gender","race",names(d)[7:length(names(d))])
# count sample size of each race
library(tidyverse)
d %>% group_by(race,type) %>% summarise(num=length(samples))
# convert wide to length
library(reshape2)
d_con <- melt(d, id.vars=c("sample","samples","type","age","gender","race"),
              variable.name="name", value.name="fpkm_log")
head(d_con)
# calculate p value
gene <- names(d)[7:length(names(d))]
p <- data.frame(matrix(c(1:(length(gene)*4)),ncol=4))
head(p)
names(p) <- c("name","pfvsm","FCfvsm","pnorvstum")
for (i in 1:length(gene)){
  female <- d_con$fpkm_log[d_con$gender=="female" & d_con$name==gene[i]]
  male <- d_con$fpkm_log[d_con$gender=="male" & d_con$name==gene[i]]
  p$pfvsm[i] <- wilcox.test(female,male)$p.value
  p$FCfvsm[i] <- mean(female,na.rm=T)-mean(male,na.rm=T)
  p$pnorvstum[i] <- wilcox.test(d_con$fpkm_log[d_con$type=="Solid Tissue Normal" & d_con$name==gene[i]],
                                d_con$fpkm_log[d_con$type=="Primary Tumor" & d_con$name==gene[i]])$p.value
}
p$name <- gene
head(p)

d_con[d_con$name==gene[i],-c(1,2,4,5,6)]
unique(d_con$name)
# plot
library("grid")
rank <- 0
pdf("20200220_luad-sex_plot3.pdf",useDingbats = F)
for(i in 1:length(gene))
{
  if(rank%%16==0){
    grid.newpage()
    rank <- 0
  }
  rank <- rank +1
  print(
    ggplot()+geom_jitter(data=na.omit(
      subset(d_con,name==gene[i])[,c("gender","fpkm_log")]),
      aes(gender,fpkm_log,colour=gender),size=0.2,height = 0)+mytheme+ggtitle(gene[i])+
      annotate("text", x=1.5,y=Inf,
               label=signif(p$pfvsm[p$name==gene[i]],digits=2), vjust=3, size=3)+
      ylab("log2(FPKM-UQ)")+scale_colour_manual(values=c("#FF7F00", "#104E8B"))+
      xlab("Gender")+ scale_x_discrete(limits=c("female","male"),label=c("Female","Male")),
    vp=viewport(.19,.21,x=(rank %% 4)*0.25+0.12,y=1-(rank %/% 4.1+0.5)*0.25)
  )
}
dev.off()

