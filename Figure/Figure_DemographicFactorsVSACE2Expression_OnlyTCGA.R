library(tidyverse)
mytheme <-  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(colour="black", size=8))+
  theme(axis.text.y = element_text(colour="black", size=8))+
  theme(axis.title=element_text(size=8))+
  #theme(plot.title=element_text(size=1.5, lineheight=.9))+
  theme(axis.ticks=element_line(colour="black",size=0.25),axis.ticks.length=unit(.2,"lines"))+
  theme(legend.title=element_text(size=8))+
  theme(legend.text=element_text(size=8))
setwd("~/projects/ACE2/FPKM-UQ_processTCGA2")
all <- read.table("Data_ACE2ExpressionInTCGAsamplesClinical.txt")
head(all)
## plot age vs ACE2 expression
head(all)
str(all)
library("lmodel2")
lm <- lmodel2(log2(subset(all,project=="LUAD" & age_at_index !="--")$FPKM+1)~
                as.numeric(as.character(subset(all,project=="LUAD" & age_at_index !="--")$age_at_index)))
cor <- cor.test(log2(subset(all,project=="LUAD" & age_at_index !="--")$FPKM+1),
                as.numeric(as.character(subset(all,project=="LUAD" & age_at_index !="--")$age_at_index)))
p1 <-ggplot(subset(all,project=="LUAD" & age_at_index !="--"),
            aes(x=as.numeric(as.character(age_at_index)),y=log2(FPKM+1)))+
  geom_point(size=0.2,colour="#00688B")+mytheme+
  xlab("Age")+
  geom_abline(intercept=as.numeric(lm$regression.results[2,2]),
              slope=as.numeric(lm$regression.results[2,3]),color='#00688B')+
  annotate("text",size=2.5, x=-Inf, y=Inf,hjust=-.4, vjust=2,label=paste("R = ",signif(cor$estimate, digit=2),sep=""))+
  annotate("text", size=2.5,x=-Inf, y=Inf,hjust=-.4, vjust=4,label=paste("P = ",signif(cor$p.value,digit=2),sep=""))
p1  

## plot gender vs ACE2 expression
gender <- as.data.frame(all %>% group_by(project,gender)%>% summarise(num=length(submitter_id)))
p2 <- ggplot(subset(all,project=="LUAD"),aes(x=gender,y=log2(FPKM+1),colour=gender))+geom_jitter(size=0.2)+
  mytheme+xlab("Gender")+ylim(-0.1,25)+
  annotate("text", x=1, y=-Inf, label=gender$num[gender$project=="LUAD"][1], vjust=-2)+
  annotate("text", x=2, y=-Inf, label=gender$num[gender$project=="LUAD"][2],  vjust=-2)+
  annotate("text", x=1.5, y=Inf, 
           label=signif(wilcox.test(log2(subset(all,project=="LUAD" & gender == "female")$FPKM+1),
                                    log2(subset(all,project=="LUAD" & gender == "male")$FPKM+1))$p.value,
                        digit=2), vjust=2)+
  scale_x_discrete(limits=c("female","male"),labels=c("Female","Male"))

p2

## plot race vs ACE2 expression
race <- as.data.frame(all %>% group_by(project,race) %>% summarise(num=length(submitter_id)))
as.data.frame(all %>% group_by(project,race2) %>% summarise(num=length(submitter_id)))
p3 <- ggplot(subset(all,project=="LUAD" & race != "american indian or alaska native" & race != "not reported"),aes(x=race,y=log2(FPKM+1),colour=race))+
  geom_jitter(size=0.2)+xlab("Race")+ylim(-0.1,25)+
  mytheme+ scale_x_discrete(limits=c("asian","black or african american","white"),labels=c("Asian", "African", "Caucasian"))+
  annotate("text", x=1, y=-Inf, label=race$num[race$project=="LUAD"][2], vjust=-2)+
  annotate("text", x=2, y=-Inf, label=race$num[race$project=="LUAD"][3],  vjust=-2)+
  annotate("text", x=3, y=-Inf, label=race$num[race$project=="LUAD"][5],  vjust=-2)+
  #annotate("text", x=4, y=-Inf, label=race$num[race$project=="LUAD"][4],  vjust=-2)+
  annotate("text", x=1.5, y=Inf, 
           label=signif(wilcox.test(log2(subset(all,project=="LUAD" & race == "asian")$FPKM+1),
                                    log2(subset(all,project=="LUAD" & race == "black or african american")$FPKM+1))$p.value,
                        digit=2), vjust=2)+
  annotate("text", x=2, y=Inf, 
           label=signif(wilcox.test(log2(subset(all,project=="LUAD" & race == "asian")$FPKM+1),
                                    log2(subset(all,project=="LUAD" & race == "white")$FPKM+1))$p.value,
                        digit=2), vjust=2)
p3
all$race2[all$race== "asian"] <- "asian"
all$race2[all$race=="black or african american"|all$race=="white"] <- "Africa and Caucasian"
wilcox.test(log2(all$FPKM[all$project=="LUAD" & all$race2=="asian"]+1),
            log2(all$FPKM[all$project=="LUAD" & all$race2=="Africa and Caucasian"] +1))

as <- mean(log2(subset(all,project=="LUAD" & race == "asian")$FPKM+1))
bl <- mean(log2(subset(all,project=="LUAD" & race == "black or african american")$FPKM+1))
wh <- mean(log2(subset(all,project=="LUAD" & race == "white")$FPKM+1))
2^(as-bl)

2^(as-wh)

library("grid")
pdf("20200202_Fig1andS1EF_jitter_only4.pdf",width=20,height=20,useDingbats = F)
print(p1,vp=viewport(.12,.095,x=.2,y=.9))
print(p2,vp=viewport(.12,.095,x=.35,y=.9))
print(p3,vp=viewport(.12,.095,x=.5,y=.9))
print(ps1,vp=viewport(.12,.095,x=.2,y=.7))
print(ps2,vp=viewport(.12,.095,x=.33,y=.7))
dev.off()
