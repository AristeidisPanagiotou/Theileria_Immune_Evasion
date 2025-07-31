setwd("/shared/projects/theileria_host_pathogen/R/")##
library(ggplot2)
library(reshape2)
library(scales)
library(ggpubr)
library(ggrepel)
library(cowplot)

##H3K18ac
chipmat <- read.csv("H3K18ac_chipmat.tsv",header = F,sep = "\t",row.names = 1)
colnames(chipmat)<-1:1200

#chipmat <- log10(chipmat+1)

#tsv <- read.csv("/shared/projects/theileria_host_pathogen/results/Kevin/mapping_hisat2/counting_featureCounts/DEA_DESeq2/Tables/dea_BL3_T_TBL3_T.tsv",sep = "\t",row.names = 1)
#tsv <- tsv[!(row.names(tsv) %in% readLines("TAgenes.txt")),]

BOLA <- grep("LOC100848815|BOLA-D",row.names(chipmat))
TLR <- grep("^TLR",row.names(chipmat))
GBP <- grep("^GBP|LOC512486$|LOC507055$|LOC511531$|LOC513659$|LOC112445995$|LOC781710",row.names(chipmat))

toplot <- rbind(melt(data.frame(chipmat[,1:600],Condition="Uninfected",Gene = row.names(chipmat))[c(BOLA,TLR,GBP),]),
                melt(data.frame(chipmat[,601:1200],Condition="Infected",Gene = row.names(chipmat))[c(BOLA,TLR,GBP),]))
#toplot$variable <- factor(as.character(gsub("^X","",toplot$variable)),levels = 1:1200)
lab <- c(rep("",19),"-3kb",rep("",279),"TSS",rep("",280),"+3kb",rep("",19))
toplot$Condition <- factor(toplot$Condition,levels = c("Uninfected","Infected"))

locs <- data.frame( LOCs = c("LOC512486","LOC781710","LOC513659","GBP5","LOC507055","GBP4","LOC112445995","LOC511531","GBP6","GBP6_1",row.names(chipmat)[c(BOLA,TLR)]),
                    Name = c("GBP1","GBP2","GBP6","GBP5","GBP4_1","GBP4","GBP7-like","GBP1_1",NA,NA,row.names(chipmat)[c(BOLA,TLR)]))
toplot$Gene <- locs[match(toplot$Gene,locs[,1]),2]
toplot <- na.omit(toplot)

toplot$group <- "Acquired Immunity"
toplot$group[grep("TLR",toplot$Gene)] <- "Innate Immunity"
toplot$group[grep("GBP",toplot$Gene)] <- "Inflammasome activation"
toplot$group <- factor(toplot$group,levels = c("Innate Immunity","Inflammasome activation","Acquired Immunity"))

ggplot(toplot, aes(x = variable, y = reorder(Gene,value), fill = value)) +
  geom_tile(width=2)+theme_bw()+xlab("")+
  theme(panel.grid = element_blank(),axis.ticks.x = element_blank(),
        strip.text = element_text(size=15), axis.text.y = element_text(hjust = 0),
        plot.title = element_text(hjust=0.5,size=15), plot.subtitle = element_text(hjust=0.5),
        aspect.ratio = 0.618,
        panel.background = element_blank())+
  viridis::scale_fill_viridis(option = "inferno",direction = -1)+
  facet_wrap("group~Condition",scales="free",ncol = 2)+
  scale_x_discrete(expand=c(-0.618,0),labels=lab)+
  labs(fill = "Norm.\nCoverage",y="")+
  ggtitle("H3K18ac intensity over promoter regions","±3kb flank")

###### H3K4me3 ######

chipmat <- read.csv("H3K4me3_chipmat.tsv",header = F,sep = "\t",row.names = 1)
colnames(chipmat)<-1:1200

#chipmat <- log10(chipmat+1)

#tsv <- read.csv("/shared/projects/theileria_host_pathogen/results/Kevin/mapping_hisat2/counting_featureCounts/DEA_DESeq2/Tables/dea_BL3_T_TBL3_T.tsv",sep = "\t",row.names = 1)
#tsv <- tsv[!(row.names(tsv) %in% readLines("TAgenes.txt")),]

BOLA <- grep("LOC100848815|BOLA-D",row.names(chipmat))
TLR <- grep("^TLR",row.names(chipmat))
GBP <- grep("^GBP|LOC512486$|LOC507055$|LOC511531$|LOC513659$|LOC112445995$|LOC781710",row.names(chipmat))

toplot <- rbind(melt(data.frame(chipmat[,1:600],Condition="Uninfected",Gene = row.names(chipmat))[c(BOLA,TLR,GBP),]),
                melt(data.frame(chipmat[,601:1200],Condition="Infected",Gene = row.names(chipmat))[c(BOLA,TLR,GBP),]))
#toplot$variable <- factor(as.character(gsub("^X","",toplot$variable)),levels = 1:1200)
lab <- c(rep("",19),"-3kb",rep("",279),"TSS",rep("",280),"+3kb",rep("",19))
toplot$Condition <- factor(toplot$Condition,levels = c("Uninfected","Infected"))

locs <- data.frame( LOCs = c("LOC512486","LOC781710","LOC513659","GBP5","LOC507055","GBP4","LOC112445995","LOC511531","GBP6","GBP6_1",row.names(chipmat)[c(BOLA,TLR)]),
                    Name = c("GBP1","GBP2","GBP6","GBP5","GBP4_1","GBP4","GBP7-like","GBP1_1",NA,NA,row.names(chipmat)[c(BOLA,TLR)]))
toplot$Gene <- locs[match(toplot$Gene,locs[,1]),2]
toplot <- na.omit(toplot)

toplot$group <- "Acquired Immunity"
toplot$group[grep("TLR",toplot$Gene)] <- "Innate Immunity"
toplot$group[grep("GBP",toplot$Gene)] <- "Inflammasome activation"
toplot$group <- factor(toplot$group,levels = c("Innate Immunity","Inflammasome activation","Acquired Immunity"))

ggplot(toplot, aes(x = variable, y = reorder(Gene,value), fill = value)) +
  geom_tile(width=2)+theme_bw()+xlab("")+
  theme(panel.grid = element_blank(),axis.ticks.x = element_blank(),
        strip.text = element_text(size=15), axis.text.y = element_text(hjust = 0),
        plot.title = element_text(hjust=0.5,size=15), plot.subtitle = element_text(hjust=0.5),
        aspect.ratio = 0.618,
        panel.background = element_blank())+
  viridis::scale_fill_viridis(option = "inferno",direction = -1)+
  facet_wrap("group~Condition",scales="free",ncol = 2)+
  scale_x_discrete(expand=c(-0.618,0),labels=lab)+
  labs(fill = "Norm.\nCoverage",y="")+
  ggtitle("H3K4me3 intensity over promoter regions","±3kb flank")

##############################
####Suppfig - Chipseq dea ####
setwd("/shared/projects/theileria_host_pathogen/R/")
lapply(c('DESeq2','ggplot2'), library, character.only=TRUE)

### statistical analysis (from featurecounts output)###

#load counts
files<-list.files("/shared/projects/theileria_host_pathogen/Deptols-TSS",".count",full.names = T)
names(files) <- gsub(".sort.bam.TSS3kb.counts","", #remove huge extension
                     list.files("/shared/projects/theileria_host_pathogen/Deptols-TSS/",".count")
)

#files<-files[grep("H3K4",names(files))][1:4] #keep duplicate H3K4me3 files
files<-files[grep("H3K18_ac",names(files))]#[c(1:2,4:5)]

tsv <- Reduce(cbind,lapply(files,function(x){
  read.table(x,header = T, sep = "\t", row.names = 1)
}))
#files are in order of previously provided gtf (3kb flanks of TSS)
colnames(tsv)<- names(files)
tsv<- tsv[!(row.names(tsv) %in% readLines("TAgenes.txt")),]


coldata <- data.frame(condition= c("Uninfected","Uninfected","Infected","Infected"),
                      row.names = colnames(tsv))
dds <- DESeqDataSetFromMatrix(countData = tsv,
                              colData = coldata,
                              design = ~ condition)

dds <- DESeq(dds)
dea <- results(dds, contrast=c("condition","Infected","Uninfected"))
dea

dea.sort<-as.data.frame(dea[order(dea$padj),])
dea.sort[,3:6]<-format(dea.sort[,3:6], digits = 5)
dea.sort[,1:2]<-round(dea.sort[,1:2],digits = 2)

dea.sort<-na.omit(dea.sort)
dea.sort[,1:6]<- sapply(dea.sort, as.numeric,USE.NAMES = T)
dea.sort<-na.omit(dea.sort)
row.names(dea.sort) <- gsub("BLA-","BOLA-",row.names(dea.sort))
#row.names(dea.sort) <- toupper(row.names(dea.sort))

write.table(dea.sort,file = "dea_H3K18_TBL3_BL3.tsv",sep = "\t",quote = F,col.names = NA)

#volcano
toplot <- as.data.frame(na.omit(dea))
tlrs <- grep("^TLR",row.names(toplot),value = T)
bolas <- grep("^BOLA-D",row.names(toplot),value = T)
gbps <- grep("^GBP|LOC512486$|LOC507055$|LOC511531$|LOC513659$|LOC112445995$|LOC781710",row.names(toplot),value = T)
keep <- c(tlrs,bolas,gbps)

toplot$Change <- "Other"
toplot[c(tlrs,gbps),"Change"] <- "TLR / GBP"
toplot[bolas,"Change"] <- "MHC Class II"

toplot$Change <- factor(toplot$Change, levels = c("TLR / GBP","MHC Class II","Other"))

#toplot <- toplot[order(-toplot$padj),]
#toplot <- rbind(toplot[grep("^TLR",row.names(toplot),invert = T),], toplot[grep("^TLR",row.names(toplot)),] ) # TLR first

toplot$txt <- NA
#keep <- row.names(toplot)[grepl("^TLR",row.names(toplot)) & abs(toplot$log2FoldChange) >1]
toplot[keep,"txt"] <- keep

toplot$size <- 2
toplot[keep,"size"] <- 3
toplot$shape <- 1
toplot[keep,"shape"] <- 19

toplot <- toplot[order(toplot$size),]

p1 <- ggplot(as.data.frame(toplot),aes(log2FoldChange,-log10(padj),color=Change,shape=shape))+scale_shape_identity()+
  geom_point(size=toplot$size)+theme_bw()+
  scale_color_manual(values=c("#000000","#A000F0","#a0a0a0"))+
  theme(panel.grid = element_blank(),axis.title = element_text(size=15),aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5),
        axis.text = element_text(size = 13), legend.position = "bottom",
        legend.text = element_text(size = 13),legend.title = element_blank() )+
  geom_vline(xintercept = c(-1,1),linetype="dashed")+geom_hline(yintercept = -log10(0.05),linetype="dashed")+
  scale_x_continuous(limits = c(-5,5),breaks = 1*(-5:5),oob = scales::oob_squish)+ #axis.text = element_text(size = 13),
  scale_y_continuous(limits = c(0,15),oob = scales::oob_squish)+
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1) ) )+
  ggtitle("ChIPSeq, Infected vs Uninfected","H3K18ac signal over gene promoters")#+  ggrepel::geom_text_repel(aes(label=txt),show.legend = FALSE)
print(p1)

#### now for h3k4me3 #########
#load counts
files<-list.files("/shared/projects/theileria_host_pathogen/Deptols-TSS",".count",full.names = T)
names(files) <- gsub(".sort.bam.TSS3kb.counts","", #remove huge extension
                     list.files("/shared/projects/theileria_host_pathogen/Deptols-TSS/",".count")
)

files<-files[grep("H3K4",names(files))][1:4] #keep duplicate H3K4me3 files
#files<-files[grep("H3K18_ac",names(files))]#[c(1:2,4:5)]

tsv <- Reduce(cbind,lapply(files,function(x){
  read.table(x,header = T, sep = "\t", row.names = 1)
}))
#files are in order of previously provided gtf (3kb flanks of TSS)
colnames(tsv)<- names(files)
tsv<- tsv[!(row.names(tsv) %in% readLines("TAgenes.txt")),]


coldata <- data.frame(condition= c("Uninfected","Uninfected","Infected","Infected"),
                      row.names = colnames(tsv))
dds <- DESeqDataSetFromMatrix(countData = tsv,
                              colData = coldata,
                              design = ~ condition)

dds <- DESeq(dds)
dea <- results(dds, contrast=c("condition","Infected","Uninfected"))
dea

dea.sort<-as.data.frame(dea[order(dea$padj),])
dea.sort[,3:6]<-format(dea.sort[,3:6], digits = 5)
dea.sort[,1:2]<-round(dea.sort[,1:2],digits = 2)

dea.sort<-na.omit(dea.sort)
dea.sort[,1:6]<- sapply(dea.sort, as.numeric,USE.NAMES = T)
dea.sort<-na.omit(dea.sort)
row.names(dea.sort) <- gsub("BLA-","BOLA-",row.names(dea.sort))
#row.names(dea.sort) <- toupper(row.names(dea.sort))

write.table(dea.sort,file = "dea_H3K4me3_TBL3_BL3.tsv",sep = "\t",quote = F,col.names = NA)

#volcano
toplot <- as.data.frame(na.omit(dea))
tlrs <- grep("^TLR",row.names(toplot),value = T)
bolas <- grep("^BOLA-D",row.names(toplot),value = T)
gbps <- grep("^GBP|LOC512486$|LOC507055$|LOC511531$|LOC513659$|LOC112445995$|LOC781710",row.names(toplot),value = T)
keep <- c(tlrs,bolas,gbps)

toplot$Change <- "Other"
toplot[c(tlrs,gbps),"Change"] <- "TLR / GBP"
toplot[bolas,"Change"] <- "MHC Class II"

toplot$Change <- factor(toplot$Change, levels = c("TLR / GBP","MHC Class II","Other"))

#toplot <- toplot[order(-toplot$padj),]
#toplot <- rbind(toplot[grep("^TLR",row.names(toplot),invert = T),], toplot[grep("^TLR",row.names(toplot)),] ) # TLR first

toplot$txt <- NA
#keep <- row.names(toplot)[grepl("^TLR",row.names(toplot)) & abs(toplot$log2FoldChange) >1]
toplot[keep,"txt"] <- keep

toplot$size <- 2
toplot[keep,"size"] <- 3
toplot$shape <- 1
toplot[keep,"shape"] <- 19

toplot <- toplot[order(toplot$size),]

p2 <- ggplot(as.data.frame(toplot),aes(log2FoldChange,-log10(padj),color=Change,shape=shape))+scale_shape_identity()+
  geom_point(size=toplot$size)+theme_bw()+
  scale_color_manual(values=c("#000000","#A000F0","#a0a0a0"))+
  theme(panel.grid = element_blank(),axis.title = element_text(size=15),aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5),
        axis.text = element_text(size = 13), legend.position = "bottom",
        legend.text = element_text(size = 13),legend.title = element_blank() )+
  geom_vline(xintercept = c(-1,1),linetype="dashed")+geom_hline(yintercept = -log10(0.05),linetype="dashed")+
  scale_x_continuous(limits = c(-5,5),breaks = 1*(-5:5),oob = scales::oob_squish)+ #axis.text = element_text(size = 13),
  scale_y_continuous(limits = c(0,15),oob = scales::oob_squish)+
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1) ) )+
  ggtitle("ChIPSeq, Infected vs Uninfected","H3K4me3 signal over gene promoters")#+  ggrepel::geom_text_repel(aes(label=txt),show.legend = FALSE)
print(p2)

### for ERC grant

toplot$RNASeq <- read.table("All_deas/dea_TBL3_K_BL3_K.tsv")[row.names(toplot),2]
toplot$shape <- 1 ; toplot["MMP9","shape"] <- 19
toplot$Change <- as.character("Other") ; toplot["MMP9","Change"] <- "MMP9" ; toplot$Change <- factor(toplot$Change, levels= c("MMP9","Other"))
toplot$txt <- NA ; toplot["MMP9","txt"] <- "MMP9"
toplot$size <- 1 ; toplot["MMP9","size"] <- 3

toplot <- toplot[order(toplot$size),]

p2 <- ggplot(as.data.frame(toplot),aes(log2FoldChange,-log10(padj),color=Change,shape=shape,label=txt))+scale_shape_identity()+
  geom_point(size=toplot$size,show.legend = F)+theme_bw()+
  geom_label(aes(x=2.25+log2FoldChange,y=2.5),size=5)+
  scale_color_manual(values=c("#000000","#a0a0a0"))+
  theme(panel.grid = element_blank(),axis.title = element_text(size=15),aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5),
        axis.text = element_text(size = 13), legend.position = "none",
        legend.text = element_text(size = 13),legend.title = element_blank() )+
  geom_vline(xintercept = c(-1,1),linetype="dashed")+geom_hline(yintercept = -log10(0.05),linetype="dashed")+
  scale_x_continuous(limits = c(-5,5),breaks = 1*(-5:5),oob = scales::oob_squish)+ #axis.text = element_text(size = 13),
  scale_y_continuous(limits = c(0,15),oob = scales::oob_squish)+
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1) ) )+
  ggtitle("ChIPSeq, Infected vs Uninfected","H3K4me3 signal over gene promoters")#+  ggrepel::geom_text_repel(aes(label=txt),show.legend = FALSE)
print(p2)

p3 <- ggplot(toplot,aes(log2FoldChange,RNASeq,color=Change,shape=shape,label=txt))+scale_shape_identity()+
  geom_point(size=toplot$size,show.legend = F)+theme_bw()+
  scale_color_manual(values=c("#000000","#a0a0a0"))+
  theme(panel.grid = element_blank(),axis.title = element_text(size=15),aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5),
        axis.text = element_text(size = 13), legend.position = "none",
        legend.text = element_text(size = 13),legend.title = element_blank() )+
  geom_vline(xintercept = 0,linetype="dashed")+   geom_hline(yintercept = 0,linetype="dashed")+
  scale_x_continuous(limits = c(-10,10),breaks = 2*(-5:5),oob = scales::oob_squish)+ #axis.text = element_text(size = 13),
  scale_y_continuous(limits = c(-10,10),breaks = 2*(-5:5),oob = scales::oob_squish)+ #axis.text = element_text(size = 13),
  xlab("H3K4me3 Log2FC")+ ylab("RNASeq Log2FC")+
  geom_label(aes(x=3+log2FoldChange, y= 1+RNASeq),size=5)+
  geom_smooth(color="green",se = F,show.legend = F)+
  annotate("text",x=7,y=-9,size=5,
           label=paste("ρ = ",round(cor(toplot$log2FoldChange,toplot$RNASeq,use = "complete.obs",method="spearman"),3)))+
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1) ) )#+  ggrepel::geom_text_repel(aes(label=txt),show.legend = FALSE)
print(p3)

cowplot::plot_grid(p2,p3,align = "hv",axis="tblr")

##########DER JUNK############
##############################
##############################
#######
toplot <- melt(data.frame(Ctrl = rowSums(chipmat[BOLA,1:600]),
                          Infected = rowSums(chipmat[BOLA,601:1200]),
                          group = "MHC Class II"))

toplot1 <- melt(data.frame(Ctrl = rowSums(chipmat[TLR,1:600]),
                           Infected = rowSums(chipmat[TLR,601:1200]),
                           group = "Toll-like Receptors"))
toplot <- rbind(toplot,toplot1)

p1<-ggplot(toplot,aes(variable,value))+
  geom_boxplot()+ #ylim(c(0,0.0001))+
  facet_wrap("group")+
  stat_compare_means(method = "wilcox.test",label = "p.signif",label.x = 1.5, label.y = 0.0001,hjust=0.5)+
  xlab("")+ylab("Normalized Coverage")+
  theme_bw()+theme(strip.text = element_text(size = 13),axis.text = element_text(size=15,colour = "#000000"),axis.title = element_text(size=15),
                   panel.grid = element_blank(),plot.title = element_text(hjust=0.5,size = 18),plot.subtitle = element_text(hjust=0.5,size = 13))+
  ggtitle("H3K18ac Intensity",subtitle = "TSS ±3kb")+
  scale_y_continuous(labels = scientific)

##H3K4me3
chipmat <- read.csv("H3K4me3_chipmat.tsv",header = F,sep = "\t",row.names = 1)

#chipmat <- log10(chipmat+1)

#tsv <- read.csv("/shared/projects/theileria_host_pathogen/results/Kevin/mapping_hisat2/counting_featureCounts/DEA_DESeq2/Tables/dea_BL3_T_TBL3_T.tsv",sep = "\t",row.names = 1)
#tsv <- tsv[!(row.names(tsv) %in% readLines("TAgenes.txt")),]

BOLA <- grep("LOC100848815|BOLA-D",row.names(chipmat))
TLR <- grep("^TLR",row.names(chipmat))

toplot <- melt(data.frame(Ctrl = rowSums(chipmat[BOLA,1:600]),
                          Infected = rowSums(chipmat[BOLA,601:1200]),
                          group = "MHC Class II"))

toplot1 <- melt(data.frame(Ctrl = rowSums(chipmat[TLR,1:600]),
                           Infected = rowSums(chipmat[TLR,601:1200]),
                           group = "Toll-like Receptors"))
toplot <- rbind(toplot,toplot1)

p2<-ggplot(toplot,aes(variable,value))+
  geom_boxplot()+ #ylim(c(0,0.0001))+
  facet_wrap("group")+
  stat_compare_means(method = "wilcox.test",label = "p.signif",label.x = 1.5, label.y = 0.0001,hjust=0.5)+
  xlab("")+ylab("Normalized Coverage")+
  theme_bw()+theme(strip.text = element_text(size = 13),axis.text = element_text(size=15,colour = "#000000"),axis.title = element_text(size=15),
                   panel.grid = element_blank(),plot.title = element_text(hjust=0.5,size = 18),plot.subtitle = element_text(hjust=0.5,size = 13))+
  ggtitle("H3K4me3 Intensity",subtitle = "TSS ±3kb")+
  scale_y_continuous(labels = scientific)


ggarrange(plotlist = align_plots(p1, p2, axis = "bt", align = "vh"),ncol = 2)
