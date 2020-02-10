setwd("~/projects/ACE2/Singapore/download/GIS031/")
# input tumor data
tumor <- read.table("Data_GSK_RSEM_rerun_expCounts_172_Tumor.tsv",header=T)
head(tumor[,c(1:5)])
## input gene length data
gene <- read.table("~/projects/ACE2/gene_anotation/longest_isoform_length.txt")
head(gene)
names(gene) <- c("gene","name","length")
head(gene)
# combine infomation of gene length and reads count in tumor camples
tumor <- merge(gene,tumor,by.x="name",by.y="Gene.symbol")
head(tumor[,c(1:10)])
#retrieve protein-coding gene
tumor_sub <- subset(tumor,Gene.type=="protein-coding")
head(tumor_sub[,c(1:10)])
# conver data frame from wide to long
library(reshape2)
tu_con <- melt(tumor_sub[,c(1,3,6:length(names(tumor_sub)))], id.vars=c("name","length"),
               variable.name="patient", value.name="count")
head(tu_con)
tu_con$patient <- gsub("$","_t",tu_con$patient)
# mark sample type
tu_con$type <- "Primary Tumor"
write.table(tu_con,"Data_GSK_RSEM_rerun_expCounts_172_Tumor_long.tsv")

# input normal data
normal <- read.table("Data_GSK_RSEM_rerun_expCounts_88_Normal.tsv",header=T)
head(normal[,c(1:5)])
#combine infomation of gene length and reads count in normal samples
normal <- merge(gene,normal,by.x="name",by.y="Gene.symbol")
head(normal[,c(1:10)])
## retrieve protein-coding gene
normal_sub <- subset(normal,Gene.type=="protein-coding")
head(normal_sub[,c(1:10)])
# conver data frame from wide to long
no_con <- melt(normal_sub[,c(1,3,6:length(names(normal_sub)))], id.vars=c("name","length"),
               variable.name="patient", value.name="count")
head(no_con)
no_con$patient <- gsub("$","_n",no_con$patient)
# mark sample type
no_con$type <- "Solid Tissue Normal"
write.table(no_con,"Data_GSK_RSEM_rerun_expCounts_88_Normal_long.tsv")


# input wide-to-long data
tumor <- read.table("Data_GSK_RSEM_rerun_expCounts_172_Tumor_long.tsv")
normal <- read.table("Data_GSK_RSEM_rerun_expCounts_88_Normal_long.tsv")
head(tumor)
# combine tumor and normal sample 
all <- rbind(tumor,normal)
head(all)
# plot reads count distribution
p1 <- ggplot(data=all,aes(x=log2(count+1),colour=patient))+geom_density()+mytheme_violin+
  facet_wrap(. ~ type)+ geom_vline(xintercept = 10)+geom_vline(xintercept = 9)+geom_vline(xintercept = 8)
p1
# calculate FPKM-UQ
case <- unique(all$patient)
uq_count <- vector(length=length(case))
length(case)
for(i in 1:length(case)){
  rank <- length(all$count[all$patient==case[i]])*0.75
  uq_count[i] <- sort(all$count[all$patient==case[i]])[rank]
  all$FPKM_UQ[all$patient==case[i]] <- (all$count[all$patient==case[i]] * 10^9)/(all$length[all$patient==case[i]] * uq_count[i])
}
head(all)
# FPKM-UQ distribution
p2 <- ggplot(data=all,aes(x=log2(FPKM_UQ+1),colour=patient))+geom_density()+mytheme_violin#+
  #facet_wrap(. ~ )#+ geom_vline(xintercept = 10)+geom_vline(xintercept = 9)+geom_vline(xintercept = 8)
p2
# write FPKM-UQ data
write.table(all,"Data_GSK_RSEM_rerun_expCounts_Normal_Tumor_long_FPKMUQ.tsv")



# input clinical data
#  mark tumor sample
clinicalt <- read.table("Data_GIS031.clinical.patient.txt",sep="\t",header=T)
head(clinicalt)
clinicalt$patient <- gsub("$","_t",clinicalt$PATIENT_ID)
# mark normal sample
clinicaln <- read.table("Data_GIS031.clinical.patient.txt",sep="\t",header=T)
head(clinicaln)
clinicaln$patient <- gsub("$","_n",clinicaln$PATIENT_ID)
# combine normal and tumor clinical data
clinical <- rbind(clinicalt,clinicaln)
head(clinical)
# merge FPKM-UQ and clincal data
all_c <- merge(all,clinical,by="patient")
head(all_c)
# check the sample size of tumor and normal samples
all_c[all_c$name=="A1BG",] %>% group_by(type) %>% summarise(num=length(patient))
# mark the race--Asian
all_c$race <- "asian"
# Convert gender information from uppercase to lowercase
all_c$gender<- "male"
all_c$gender[all_c$GENDER=="Female"] <- "female"
write.table(all_c,"Data_GSK_RSEM_rerun_expCounts_Normal_Tumor_long_FPKMUQ_clinical.tsv")
