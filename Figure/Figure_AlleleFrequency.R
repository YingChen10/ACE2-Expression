library(tidyverse)
####################

## 1KGP ##

####################

#Read in variations extracted from the 1000 Genome Project 
SNP<-read.table("/Dell/Dell1/shankj/projects/nCov/Extract_Variations.txt",header=T,sep = "\t")
head(SNP)
#Read in variation annotations extracted from Ensembl
Esembl<-read.csv("/Dell/Dell1/shankj/projects/nCov/Variation_annotations.csv",header=T)
nrow(Esembl)

#variations from the 1000 Genome Project to be annotated
Merge<-merge(SNP,Esembl,by.x ="variation", by.y= "Variant.ID", all.x = T)
#remove replications in non coding transcript 
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


p <- ggplot()+geom_point(data = test,aes(x = log10(AF), y = log10(EAS)), size=0.6,color="#104E8B")+
  geom_abline(intercept=0, slope=1)+
  geom_point(data = subset(test,hl=="Mis"),aes(x = log10(AF), y = log10(EAS)), size=0.6,color="#FF4040")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  ylab("log10(Allele frequency of EAS populations)")+
  xlab("log10(Allele frequency of all populations)")+  
  scale_x_continuous(breaks = c(-4,-3,-2,-1,0),limits = c(-4,0))+  
  scale_y_continuous(breaks = c(-4,-3,-2,-1,0),limits = c(-4,0))
p

#correlation test
cor.test(test$AF,test$EAS,method = "pearson")

#Save
library("grid")
pdf("fig1-3.pdf",width=20,height=20,useDingbats = F)
print(p,vp=viewport(.1,.1,x=.2,y=.9))
dev.off()


####################

## gnomAD ##

####################

#Read in variations extracted from the gnomAD 
gnomAD<-read_csv("/Dell/Dell1/shankj/projects/nCov/gnomAD_v3_ENSG00000130234_2020_02_20_10_54_05.csv")
head(gnomAD)
colnames(gnomAD)
gnomAD$Pos<-paste0("X:",gnomAD$Position)
gnomAD$new <- gnomAD$`Allele Count East Asian`

#if there is no alternate allele count in East Asian, then the count of alternate alleles is assigned to 0.5
#it is a convention
gnomAD$new[gnomAD$`Allele Count East Asian` == 0] <- 0.5

#allele frequence
gnomAD$Freq_EastAsian<-gnomAD$new/gnomAD$`Allele Number East Asian`

#Select only East Asian
gnomAD_Plot<-select(gnomAD, Pos,rsID,Reference,Alternate,`Allele Frequency`,Freq_EastAsian,Annotation)
colnames(gnomAD_Plot)<-c("Pos","rsID","Reference","Alternate","All","East_Asian","Type")
head(gnomAD_Plot)
unique(gnomAD_Plot$Type)

#Plot 
gnomAD_Plot$hl <- "Others"
gnomAD_Plot$hl[gnomAD_Plot$Type %in% 
                 c("missense_variant",
                   "frameshift_variant")] <- "Mis"
unique(gnomAD_Plot$hl)
gnomAD_Plot %>% group_by(hl) %>% summarise(num=length(Pos))
length(gnomAD_Plot$hl[gnomAD_Plot$hl=="Mis"])

p <- ggplot()+geom_point(data = gnomAD_Plot,aes(x = log10(All), y = log10(East_Asian)), size=0.6,color="#104E8B")+
  geom_abline(intercept=0, slope=1)+
  geom_point(data = subset(gnomAD_Plot,hl=="Mis"),aes(x = log10(All), y = log10(East_Asian)), size=0.6,color="#FF4040")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  ylab("log10(Allele frequency of EAS populations)")+
  xlab("log10(Allele frequency of all populations)")+  
  scale_x_continuous(breaks = c(-5,-4,-3,-2,-1,0),limits = c(-5.5,0))+  
  scale_y_continuous(breaks = c(-5,-4,-3,-2,-1,0),limits = c(-5.5,0))
p

#correlation test
cor.test(gnomAD_Plot$All,gnomAD_Plot$East_Asian,method = "pearson")

pdf("gnomAD.pdf",width=20,height=20,useDingbats = F)
print(p,vp=viewport(.1,.1,x=.2,y=.9))
dev.off()