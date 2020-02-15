setwd("~/projects/ACE2/FPKM-UQ_processTCGA2")
# input 
ace2 <- read.table("20200205_combine/Data_ACE2ExpressionNormalizedByMedian.txt")
head(ace2)
# calculate SD and variation
sd(ace2$nor_median,na.rm=T)
var(ace2$nor_median)
# remark age and race
ace2$AGE <- as.numeric(as.character(ace2$AGE))
ace2$RACE[ace2$RACE == "american indian or alaska native" | ace2$RACE == "not reported"] <-  NA
# anova all race
summary(a <- lm(ace2$nor_median ~ ace2$AGE + ace2$RACE + ace2$GENDER))
anova(a)

# plot age vs ACE2 expression
# correlation between age and ACE2 expression
lm <- lmodel2(ace2$nor_median~as.numeric(as.character(ace2$AGE)))
lm
cor <- cor.test(ace2$nor_median,as.numeric(as.character(ace2$AGE)))
cor
# calculate sample size of sample with reported age
n_age <- length(na.omit(ace2[,c("nor_median","AGE")])$AGE)

p1 <-ggplot(data=na.omit(ace2[,c("nor_median","AGE")]),aes(x=AGE,y=nor_median))+
  geom_point(size=0.2,colour="#00688B")+mytheme_violin+ylim(-0.1,25)+
  xlab("Age")+
  geom_abline(intercept=as.numeric(lm$regression.results[2,2]),
              slope=as.numeric(lm$regression.results[2,3]),color='#00688B')+
  annotate("text",size=2, x=-Inf, y=Inf,hjust=-.4, vjust=2,label=paste("R = ",signif(cor$estimate, digit=2),sep=""))+
  annotate("text", size=2,x=-Inf, y=Inf,hjust=-.4, vjust=4,label=paste("P = ",signif(cor$p.value,digit=2),sep=""))+
  annotate("text", size=2, x=-Inf, y=Inf,hjust=-.4, vjust=6,label=paste("N = ",n_age,sep=""))
p1  

# plot gender vs ACE2 expression
ggplot(ace2,aes(x=GENDER,y=nor_median))+geom_boxplot()
ggplot(ace2,aes(x=GENDER,y=nor_median))+geom_boxplot()
wilcox.test(ace2$nor_median[ace2$GENDER=="male"],ace2$nor_median[ace2$GENDER=="female"])

ggplot(subset(ace2,RACE != "american indian or alaska native" & RACE != "not reported"),
       aes(x=GENDER,y=nor_median))+geom_boxplot(notch=T)
wilcox.test(subset(ace2,GENDER=="male" & RACE != "american indian or alaska native" & RACE != "not reported")$nor_median,
            subset(ace2,GENDER=="female" & RACE != "american indian or alaska native" & RACE != "not reported")$nor_median)

# calculate sample size of each gender
unique(ace2$GENDER)
gender <- as.data.frame(ace2 %>% group_by(GENDER)%>% summarise(num=length(patient)))
gender
p2 <- ggplot(data=ace2,aes(x=GENDER,y=nor_median,colour=GENDER))+geom_jitter(size=0.2)+
  mytheme_violin+xlab("Sex")+ylim(-0.1,25)+
  annotate("text", x=1, y=-Inf, label=gender$num[1], vjust=-2)+
  annotate("text", x=2, y=-Inf, label=gender$num[2],  vjust=-2)+
  annotate("text", x=1.5, y=Inf, 
           label=signif(wilcox.test(ace2$nor_median[ace2$GENDER=="male"],
                                                ace2$nor_median[ace2$GENDER=="female"])$p.value,digit=2), vjust=2)+
  scale_x_discrete(limits=c("female","male"),labels=c("Female","Male"))
p2

# calculate sample size in each race
race <- ace2 %>% group_by(RACE) %>% summarise(num=length(patient))
race
# plot race vs ACE2 expression
p3 <- ggplot(subset(ace2,RACE != "american indian or alaska native" & RACE != "not reported"),
             aes(x=RACE,y=nor_median,colour=RACE))+
  geom_jitter(size=0.2)+xlab("Race")+ylim(-0.1,25)+
  mytheme_violin+ 
  scale_x_discrete(limits=c("asian","black or african american","white"),labels=c("Asian", "African", "European"))+
  annotate("text", x=1, y=-Inf, label=race$num[1], vjust=-2)+
  annotate("text", x=2, y=-Inf, label=race$num[2],  vjust=-2)+
  annotate("text", x=3, y=-Inf, label=race$num[3],  vjust=-2)+
  annotate("text", x=1.5, y=Inf, 
           label=signif(wilcox.test(subset(ace2, RACE == "asian")$nor_median,
                                    subset(ace2, RACE == "black or african american")$nor_median)$p.value,
                        digit=2), vjust=2)+
  annotate("text", x=2, y=Inf, 
           label=signif(wilcox.test(subset(ace2, RACE == "asian")$nor_median,
                                    subset(ace2,RACE == "white")$nor_median)$p.value,
                        digit=2), vjust=2)
p3
library("grid")
pdf("20200205_combine/20200210_combine_fig1.pdf",width=20,height=20,useDingbats = F)
print(p1,vp=viewport(.1,.1,x=.2,y=.9))
print(p2,vp=viewport(.1,.1,x=.31,y=.9))
print(p3,vp=viewport(.1,.1,x=.42,y=.9))
dev.off()
