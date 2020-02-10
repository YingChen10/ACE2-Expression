setwd("~/projects/ACE2/FPKM-UQ_processTCGA2")
# input TCGA LUAD sample all gene expression
luad1 <- readRDS("ACE2_expr_cor/Data_TCGALUADSampleExpression1.rds")
luad2 <- readRDS("ACE2_expr_cor/Data_TCGALUADSampleExpression2.rds")
head(luad1)
luad_tcga <- rbind(luad1,luad2)
# remark patient id with  tumor or normal
luad_tcga_t <- subset(luad_tcga,tissue=="Primary Tumor"|tissue=="Recurrent Tumor")
luad_tcga_n <- luad_tcga[luad_tcga$tissue=="Solid Tissue Normal",]  
luad_tcga_t$case  <- gsub("$","_t",luad_tcga_t$case)
luad_tcga_n$case  <- gsub("$","_n",luad_tcga_n$case) 
# combine normal and tumor
luad_tcga_combine <- rbind(luad_tcga_t,luad_tcga_n)
##saveRDS(luad_tcga_combine,"20200205_combine/TCGALUADGeneExpression.rds")
head(luad_tcga)
# subset TCGA LUAD
luad_tcga <- luad_tcga_combine[,c("name","case","tissue","expr")]
head(luad_tcga)
names(luad_tcga) <- c("name","patient","type","FPKM_UQ")
luad_tcga$SOURCE <- "TCGA"
# input OncoSG smaple
OncoSG <- read.table("~/projects/ACE2/Singapore/download/GIS031/GSK_RSEM_rerun_expCounts_Normal_Tumor_long_FPKMUQ.tsv")
head(OncoSG)
# subset OncoSG LUAD
OncoSG <- OncoSG[,c("name","patient","type","FPKM_UQ")]
head(OncoSG)
OncoSG$SOURCE <- "OncoSG"
# combine TCGA and OncoSG
all <- rbind(OncoSG,luad_tcga)
head(all) 
median_TCGA <- as.data.frame(all %>% group_by(SOURCE) %>% summarise(median = median(log2(FPKM_UQ+1))))[2,2]
median_TCGA 
median_OncoSG <- as.data.frame(all %>% group_by(SOURCE) %>% summarise(median = median(log2(FPKM_UQ+1))))[1,2]
median_OncoSG
# normaliezed by median combine
all$nor_median <- log2(all$FPKM_UQ+1)
all$nor_median[all$SOURCE=="OncoSG"] <- log2(all$FPKM_UQ[all$SOURCE=="OncoSG"]+1)+(median_TCGA-median_OncoSG)
# write table
write.table(all,"20200205_combine/all_normalized_by_median.txt")
write.table(subset(all,name=="ACE2"),"20200205_combine/ACE2_normalized_by_median.txt")
## figS4
p2 <- ggplot(data=all,aes(x=nor_median,colour=patient))+geom_density(alpha=0.1)+mytheme_violin+
  scale_colour_manual(values=c(rep("#EE7942",260), rep("#104E8B",594)))
p2

p1 <- ggplot(data=all,aes(x=log2(FPKM_UQ+1),colour=patient))+geom_density(alpha=0.1)+mytheme_violin+
  scale_colour_manual(values=c(rep("#EE7942",260), rep("#104E8B",594)))
p1

pdf("20200205_combine/20200202_FigS4.pdf",width=20,height=20,useDingbats = F)
print(p1,vp=viewport(.12,.09,x=.2,y=.8))
print(p2,vp=viewport(.12,.09,x=.4,y=.8))
dev.off()


