
mydata=read.csv("IMGVR_diamond_1000.tsv",header = FALSE, sep = "\t")
mydata$V2 <-gsub("\\|.*","",mydata$V2)
head(mydata)
write.table(mydata,"IMGVR_diamond.rename.tsv",quote = F,sep = "\t",row.names = FALSE,col.names = FALSE)

library(dplyr)
tax_data=read.table("IMGVR_diamond.rename.tsv",header=T)

tax_dat1 <- tax_data %>%
  group_by(qseqid) %>%
  filter(bitscore==max(bitscore))
tax_dat2<- tax_dat1 %>%
  group_by(qseqid) %>%
  filter(pident==max(pident))
write.table(tax_dat2,"IMGVR_diamond.rename.filtered.tsv",quote = F,row.names = FALSE,sep = "\t")

library(stringr)
a=str_count(tax_dat2$Taxonomic.classification)
b=data.frame(a)
tax_data2_new=data.frame(tax_dat2,b)

tax_dat3 <- tax_data2_new %>%
  group_by(qseqid) %>%
  filter(a==max(a))
tax_dat3_durep=tax_dat3[!duplicated(tax_dat3$qseqid), ]

write.csv(tax_dat3_durep,"taxonomy_new_unique.csv",quote = F)

data2_n <- data2 %>%
  group_by(sseqid) %>%
  filter(bit.score==max(bit.score))