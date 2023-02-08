#PcoA
library(vegan)
library(ggplot2)
setwd("/Users/shicong/Desktop/Pangenome/Pangenome/Diversity")
otu=read.csv("ARG_RPKM_coverage.csv",header = T,row.names = 1)
otu=t(otu)
G=read.csv("group.csv",header = T, row.names = 1) 
otu.bry=vegdist(otu) 
otu.pcoa=cmdscale(otu.bry,eig=T)
point= otu.pcoa$points
colnames(point)=c("PCoA1","PCoA2")
pcoa=data.frame(point,G) 
eig=otu.pcoa$eig
format(100 * eig[1] / sum(eig))
format(100 * eig[2] / sum(eig))
ggplot(pcoa,aes(x=PCoA1,y=PCoA2,colour=group))+geom_point(size=2.5) + xlab("PCoA1 (24.63%)") + ylab("PCoA1 (17.52%)") +scale_shape_manual(values = c(8,0,1,10,3,4,7,3,11,5,2))+ scale_colour_manual(values=c("#C9716A","#EBBB46","#5FA35D","#6490C3"))+gghighlight::gghighlight(label_key = Sample_type)+facet_wrap(vars(Sample_type))

stat_ellipse(aes(fill=dm_diagnosed),type="norm",geom="polygon",alpha=0.2,color=NA)
##PERMANOVA
 
adonis(otu.bry~G$dm_diagnosed,permu=999)
##anosim
anosim(otu.bry,G$dm_diagnosed,permutations = 999)  
#rarecurve
setwd("/Users/shicong/Desktop/Pangenome/Pangenome/virus_coverage/rarefraction")

rarecurve(otu)