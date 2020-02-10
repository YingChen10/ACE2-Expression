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
## plot pathologic stage vs ACE2 expression
all$pathologic_stage_rank[all$ajcc_pathologic_stage %in% "Stage I"|all$ajcc_pathologic_stage %in% "Stage IA"|all$ajcc_pathologic_stage %in% "Stage IB"|all$ajcc_pathologic_stage %in% "Stage IC"] <- "Stage I"
all$pathologic_stage_rank[all$ajcc_pathologic_stage %in% "Stage II"|all$ajcc_pathologic_stage %in% "Stage IIA"|all$ajcc_pathologic_stage %in% "Stage IIB"|all$ajcc_pathologic_stage %in% "Stage IIC"] <- "Stage II"
all$pathologic_stage_rank[all$ajcc_pathologic_stage %in% "Stage III"|all$ajcc_pathologic_stage %in% "Stage IIIA"|all$ajcc_pathologic_stage %in% "Stage IIIB"|all$ajcc_pathologic_stage %in% "Stage IIIC"] <- "Stage III"
all$pathologic_stage_rank[all$ajcc_pathologic_stage %in% "Stage IV"] <- "Stage IV"
as.data.frame(all %>% group_by(project,ajcc_pathologic_stage) %>% summarise(num=length(submitter_id)))
stage <-as.data.frame(all %>% group_by(project,pathologic_stage_rank) %>% summarise(num=length(submitter_id)))

ps1 <- ggplot(data=subset(all,project=="LUAD" & ajcc_pathologic_stage != "--"),aes(x=pathologic_stage_rank,y=log2(FPKM+1),colour=pathologic_stage_rank))+geom_jitter(size=0.2)+
  mytheme+xlab("Pathologic stage")+
  annotate("text", x=1, y=-Inf, label=stage$num[stage$project=="LUAD"][1], vjust=-2)+
  annotate("text", x=2, y=-Inf, label=stage$num[stage$project=="LUAD"][2],  vjust=-2)+
  annotate("text", x=3, y=-Inf, label=stage$num[stage$project=="LUAD"][3], vjust=-2)+
  annotate("text", x=4, y=-Inf, label=stage$num[stage$project=="LUAD"][4],  vjust=-2)
ps1

ps2 <- ggplot(data=subset(all,project=="LUSC" & ajcc_pathologic_stage != "--"),aes(x=pathologic_stage_rank,y=log2(FPKM+1),colour=pathologic_stage_rank))+geom_jitter(size=0.2)+
  mytheme+xlab("Pathologic stage")+
  annotate("text", x=1, y=-Inf, label=stage$num[stage$project=="LUSC"][1], vjust=-2)+
  annotate("text", x=2, y=-Inf, label=stage$num[stage$project=="LUSC"][2],  vjust=-2)+
  annotate("text", x=3, y=-Inf, label=stage$num[stage$project=="LUSC"][3], vjust=-2)+
  annotate("text", x=4, y=-Inf, label=stage$num[stage$project=="LUSC"][4],  vjust=-2)
ps2
