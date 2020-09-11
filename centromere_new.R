# Read table from file. The table has 2 columns: geneID and GO number (GO:0000001)
.libPaths("D:\\Programming\\R_working\\Rpackages")
setwd("D:\\Programming\\R_working\\Rworking\\ALL_BGI")
library(topGO)
library(data.table)
library(Rgraphviz)
library(xlsx)
library(ggplot2)
library(ggsci)



annotation_AC = fread("CLUSTER_TABLE_AC.csv_parsed")
annotation_AF= fread("CLUSTER_TABLE_AF.csv_parsed")

annotation_AC[Supercluster == 1]

#fwrite(annotation_AC, file = "formated_CLUSTER_TABLE_AC.csv", sep = "\t", row.names = F, quote = F)
unique(annotation_AC$Supercluster_best_hit)

AC_reads = 2412538
annotation_AC$portion = annotation_AC$Size*100/AC_reads
annotation_AC$Cluster = paste("CL",annotation_AC$Cluster, sep = "")


AF_TR_vs_contig = fread("TR_vs_contigs_AF_shortNames")#fread("TR_vs_contigs_AF.tab")
AF_TR_vs_contig$cluster = unname(sapply(AF_TR_vs_contig$V2, FUN=function(x){strsplit(x,'Contig')[[1]][1]}))



AF_reads = 1282472
annotation_AF$portion = annotation_AF$Size*100/AF_reads
annotation_AF$Cluster = paste("CL",annotation_AF$Cluster, sep = "")

AC_TR_vs_contig = fread("TR_vs_contigs_AC_shortNames")
AC_TR_vs_contig$cluster = unname(sapply(AC_TR_vs_contig$V2, FUN=function(x){strsplit(x,'Contig')[[1]][1]}))


annotation_AF$Species = "AF"
annotation_AC$Species = "AC"


##Satellitome analysis
byCluster_AC = AC_TR_vs_contig[,.(TRs = paste(unique(V1), collapse = ","), count = length(V2)), by = "cluster"]
byCluster_AF = AF_TR_vs_contig[,.(TRs = paste(unique(V1), collapse = ","), count = length(V2)), by = "cluster"]

### add to the annotation table
annotation_AC = merge(annotation_AC, all.x = T, by.x = "Cluster", byCluster_AC, by.y = 'cluster', all.y = T)
annotation_AC[!(is.na(TRs))]$SuperFamily = "Satellite"


annotation_AF = merge(annotation_AF, all.x = T, by.x = "Cluster", byCluster_AF, by.y = 'cluster', all.y = T)
annotation_AF[!(is.na(TRs))]$SuperFamily = "Satellite"

annotation_AC[is.na(Supercluster)]
annotation_AC[TRs == "TR2CL137"]$portion =  0.014
annotation_AC[TRs == "TR2CL137"]$Species =  "AC"

annotation_AC = annotation_AC[!(is.na(portion))]

annotation_AC[!(is.na(TRs))]$Lineage = annotation_AC[!(is.na(TRs))]$TRs
annotation_AF[!(is.na(TRs))]$Lineage = annotation_AF[!(is.na(TRs))]$TRs
annotation_AF = annotation_AF[!(is.na(portion))]

##plot SuperFamily barplot

joined_annotaion = rbind(annotation_AC,annotation_AF)
joined_annotaion = as.data.table(joined_annotaion)

joined_annotaion[SuperFamily == "Class_I"]$SuperFamily = "LINE"
joined_annotaion[SuperFamily == "TIR"]$SuperFamily = "Class II TE"


##FIgure 1a
ggplot(joined_annotaion[SuperFamily != "organelle",.(portion = sum(portion)), by = c("SuperFamily", "Species")]) + 
  geom_bar(aes(SuperFamily, portion, fill = Species), position="dodge",stat='identity') +
  #ggtitle("Comparative analysis of repeatome of three Allium species") + 
  #coord_flip() + 
  ylab("Genome proportion") +
  xlab("Repeat type") +
  theme_bw() + 
  theme(axis.text = element_text(color="black", size=24),
        text = element_text(color="black",size=24),
        axis.title = element_text(face = "bold")) +
  scale_fill_aaas(name = "Species", labels = c("AC", "AF"))  #+ facet_grid(~variable, scales = 'free') 

##plot repeat leneage barplot
joined_annotaion[,.(portion = sum(portion)), by = c("Lineage", "Species")]


##FIgure 2b
ggplot(joined_annotaion[SuperFamily  != "Satellite" & 
                          Lineage  != "mitochondria" & 
                          Lineage  != "plastid" &
                        Lineage  != "Unknown", 
                        .(portion = sum(portion)), by = c("Lineage", "Species")]) + 
  geom_bar(aes(Lineage, portion, fill = Species), position="dodge",stat='identity') +
  #ggtitle("Comparative analysis of repeatome of three Allium species") + 
  #coord_flip() + 
  ylab("Genome proportion") +
  xlab("Mobile elements") +
  theme_bw() + 
  theme(axis.text.y = element_text(color="black", size=24),
        axis.text.x = element_text(color="black", size=24, angle = 45, hjust = 1),
        text = element_text(color="black",size=24),
        axis.title = element_text(face = "bold")) +
  scale_fill_aaas(name = "Species", labels = c("AC", "AF"))  #+ facet_grid(~variable, scales = 'free') 



mer_satell_AC = merge(byCluster_AC, by.x = "cluster",annotation_AC[,c(1,6)], by.y = "Cluster", all.x = T)
mer_satell_AC[TRs == "TR2CL137"]$portion =  0.014
perTR_portion_AC = mer_satell_AC[!(is.na(portion)),.(portion = sum(portion)), by = "TRs"]
perTR_portion_AC$Species = "AC"

mer_satell_AF = merge(byCluster_AF, by.x = "cluster",annotation_AF[,c(1,6)], by.y = "Cluster", all.x = T)
perTR_portion_AF = mer_satell_AF[!(is.na(portion)),.(portion = sum(portion)), by = "TRs"]
perTR_portion_AF$Species = "AF"

mer_satellite_ACAF = rbind(perTR_portion_AF,perTR_portion_AC)

#Figure 2c
ggplot(mer_satellite_ACAF) + 
  geom_bar(aes(as.character(TRs), portion, fill = Species), stat='identity') +
  #ggtitle("Comparative analysis of repeatome of three Allium species") + 
  #coord_flip() + 
  ylab("Genome proportion") +
  xlab("Tandem repeat") +
  theme_bw() + 
  theme(axis.text = element_text(color="black", size=18),
        axis.text.x = element_blank(),
        text = element_text(color="black",size=18),
        axis.title = element_text(face = "bold")
        ) + 
  scale_fill_aaas() +
  facet_wrap(.~TRs, scales = "free")#name = "Species", labels = c("AC", "AF"))  #+ 

######

