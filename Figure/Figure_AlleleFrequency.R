library(tidyverse)

#Read in variations extracted from the 1000 Genome Project 
SNP<-read.table("/Dell/Dell1/shankj/projects/nCov/Extract_Variations.txt",header=T,sep = "\t")
head(SNP)
#Read in variation annotations extracted from Ensembl
Esembl<-read.csv("/Dell/Dell1/shankj/projects/nCov/Variation_annotations.csv",header=T)
nrow(Esembl)

#variations from the 1000 Genome Project to be annotated
Merge<-merge(SNP,Esembl,by.x ="variation", by.y= "Variant.ID", all.x = T)
#to filter replicated variantions in non-coding transcripts 
test<-unique(Merge %>% filter(!grepl('non coding', Conseq..Type)) %>% select(variation, Pos, AF,EAS,Conseq..Type))
head(test)
nrow(test)

#Plot 
names(test) <- c("SNP","Pos","AF","EAS","type")
test$hl <- "Others"
test$hl[test$type %in% 
          c("missense variant",
            "missense variant~splice region variant",
            "stop gained")] <- "Mis"

unique(test$hl)
test %>% group_by(hl) %>% summarise(num=length(Pos))

p <- ggplot()+geom_point(data = test,aes(x = 100*AF, y = 100*EAS),size=0.6,color="#104E8B")+
  geom_point(data = subset(test,hl=="Mis"),aes(x = 100*AF, y = 100*EAS), size=0.6,color="#FF4040")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  ylab("Allele frequency of EAS populations (%)")+
  xlab("Allele frequency of all populations (%)")+
  geom_abline(intercept=0, slope=1)
p

library("grid")
pdf("fig1-3.pdf",width=20,height=20,useDingbats = F)
print(p,vp=viewport(.1,.1,x=.2,y=.9))
dev.off()
