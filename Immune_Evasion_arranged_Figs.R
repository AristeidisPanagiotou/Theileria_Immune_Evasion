
#library(ggpubr)
########## FIGURE 1 ##########

pt <- "C:/Users/Aris/gsea_home/output/feb02/KMGW_B_cancer_inf_vs_un.Gsea.1706871556397"
#grab all gsea reports from treated vs untreated and
files<-  list.files(pt,pattern = "^gsea_report.+tsv", recursive = T,full.names = T)

gseas<-rbind(read.csv(files[1],sep="\t")[,c("NAME","NES","FDR.q.val")],
             read.csv(files[2],sep="\t")[,c("NAME","NES","FDR.q.val")])

#preprocess for plotting (x shape for non significant, log the qval, shorten names)
gseas<-gseas[!gseas$NES=="---",]
gseas$NES <- as.numeric(gseas$NES)
gseas$shape<- 1
gseas[gseas$FDR.q.val<0.05,"shape"] <- 19
gseas$`-log10(qval)` <- -log10(gseas$FDR.q.val + 1e-3)

sapply(gseas, class) #check classes
gseas$xNAME <- gsub("_"," ",gseas$NAME)
gseas$xNAME <- paste0(toupper(substr(gseas$xNAME, 1, 1)), tolower(substr(gseas$xNAME, 2, nchar(gseas$xNAME))))


ggplot(gseas, aes(1,reorder(xNAME,`-log10(qval)`)))+ scale_shape_identity()+
  geom_point(aes(size = `-log10(qval)`,colour= NES,shape= shape), show.legend = TRUE)+
  theme_bw()+
  ggtitle("RNASeq GSEA summary","Cancer Hallmarks in infected cells") + 
  theme(plot.title = element_text(hjust=0.5,size = 18),
        plot.subtitle = element_text(hjust=0.5,size=15),
        legend.title = element_text(size=12),legend.text = element_text(size=10),
        panel.grid = element_blank(),
        aspect.ratio = 4,axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12,vjust = 0.3),
        axis.ticks.x=element_blank(),axis.title = element_blank()) +
  scale_color_gradient2(low="blue", mid="white", high="red", 
                        limits = c(-2, 2), oob = scales::squish)+
  scale_size_continuous(breaks = 0:3,limits = c(0, 3),range = c(0, 10) )#+  geom_hline(yintercept = c(8.5,6.5),linetype="longdash")


toload <- grep(paste(gseas$NAME,collapse="|"),list.files(pt,"\\.tsv$",full.names = T),value = T)
names(toload) <- grep(paste(gseas$NAME,collapse="|"),list.files(pt,"\\.tsv$"),value = T)

tmp <- sapply(names(toload),function(x){
names(toload)[1]
  j <- read.csv(toload[x],sep="\t")[,c(2,5,7)]
  #j <- j[c(1:20,(nrow(j)-19):nrow(j)),]
  j <- j[c(1:10,(nrow(j)-9):nrow(j)),]
  j <- j[j$CORE.ENRICHMENT == "Yes",]
  j <- cbind(j,gsub("\\.tsv$","",x))
  j
},simplify = F)
names(tmp) <- gsub("\\.tsv$","",names(tmp))

tmp <- tmp[c(8,7,1,4)]
tmp <- Reduce(rbind,tmp)
colnames(tmp)[4] <- "Geneset"

###AUUUUUUUUUGH retranslate genesymbols to bos ##########
chip <- read.csv("C:/Users/Aris/Desktop/GSEA/Inputs/BT_HS.chip",sep="\t")[,-3] ; row.names(chip)<-chip[,2]
cbind(tmp$SYMBOL,toupper(chip[tmp$SYMBOL,1])) #check
tmp$SYMBOL<-toupper(chip[tmp$SYMBOL,1])


tsv <- read.csv("C:/Users/Aris/Desktop/GSEA/Inputs/KMGW_B_inf_un.gct", sep="\t",skip = 2)[,-2]
tsv <- tsv[tsv$NAME %in% tmp$SYMBOL,]
row.names(tsv) <- tsv$NAME

toplot <- reshape::melt(data.frame(Gene=tmp$SYMBOL,tsv[tmp$SYMBOL,-1],Geneset=tmp$Geneset))


colnames(toplot)[3:4]<-c("Sample","TPM")

cp <- coord_polar(theta = "x"); cp$is_free <- function() TRUE

toplot$Sample <- factor(toplot$Sample,levels = rev(unique(toplot$Sample)))
toplot$Geneset <- gsub("_"," ",toplot$Geneset)
toplot$Geneset <- paste0(toupper(substr(toplot$Geneset, 1, 1)), tolower(substr(toplot$Geneset, 2, nchar(toplot$Geneset))))

##calculate first 10
#lvs <- unlist(lapply(unique(toplot$Geneset),function(x){
#  tmpz <- toplot[toplot$Geneset == x,]
#  levels(reorder(tmpz$Gene,-tmpz$TPM))[1:10]
#}))
#toplot <- toplot[toplot$Gene %in% lvs,]

ggplot(toplot, aes(Sample,reorder(Gene,TPM),fill=log10(TPM+1)))+
  geom_tile(show.legend = TRUE)+
  theme_bw()+ylab("Gene")+
  viridis::scale_fill_viridis(option = "inferno",direction = -1)+ #inferno/mako
  #scale_fill_gradient2(low="blue", mid="white", high="red", midpoint = 2,limits = c(0, 4), oob = scales::squish)+
  ggtitle("Top 10 Core Enrichment genes") + 
  theme(plot.title = element_text(hjust=0.5,size = 18),panel.grid = element_blank(),
        aspect.ratio = 1,legend.position = "right",panel.border = element_blank(),
        strip.text = element_text(size = 12),
        axis.ticks.x=element_blank(),
        axis.ticks.length.y = unit(-3.6,"cm"),strip.background = element_blank(),
        axis.ticks.y = element_line(linetype = "dotted",linewidth = 0.5,colour="#ffffff"), #,colour="#D0D0D0"
        axis.title = element_blank(),
        axis.text.x = element_text(face = "bold",size = 15,colour = "#000000"),
        axis.text.y = element_text(size = 9,colour = "#000000"))+
  scale_x_discrete(labels=c(rep("",12),"\nInfected    ",rep("",4),"\nUninfected ",rep("",12)))+
  cp+
  facet_wrap("Geneset",scales = "free",ncol = 2,nrow = 2)
#save

dats <- c("Villares et al. TBL3","Cheeseman  et al. TBL3","Rchiad et al. TBL3","Rchiad et al. TBL20","KW Unpublished  TBL20",
          "Villares et al. BL3" ,"Cheeseman  et al. BL3" ,"Rchiad et al. BL3" ,"Rchiad et al. BL20" ,"KW Unpublished  BL20")



colorz <- c(rep("#FF0000",3),rep("#FF4040",3),rep("#FF6060",3),rep("#FF8080",3),rep("#FFA0A0",3),
            rep("#0000FF",3),rep("#4040FF",3),rep("#6060FF",3),rep("#8080FF",3),rep("#A0A0FF",3))


lb <- data.frame(lbs=dats,x=3*0:9 + 2,y=22)


ggplot(toplot, aes(Sample,Gene,fill=Sample,color=Sample))+
  geom_tile(height=1.2,width = 1.05, show.legend = F)+
  theme_bw()+ylab("Gene")+
  scale_fill_manual(values=colorz)+
  scale_color_manual(values=colorz)+
  #ggtitle("Top 20 Core Enrichment genes") + 
  theme(plot.title = element_text(hjust=0.5),plot.subtitle = element_text(hjust=0.5),panel.grid = element_blank(),
        aspect.ratio = 1,strip.text = element_text(size = 12),
        axis.ticks=element_blank(),axis.title = element_blank(),
        panel.border = element_blank(),
        axis.text= element_blank())+
  coord_polar(theta = "x") +
  annotate("text",x=3*0:9 + 2,22,label=dats,angle=c(432-36*0:4,-432+36*4:0),size = 4.5,fontface="bold")+

  geom_vline(xintercept = 0.5+3*0:9,linewidth=1.02)+
  geom_hline(yintercept = 38.5,linewidth=1.1)
 
########## FIGURE 2 ##########
########## Volcano  ##########

dea <- read.csv("dea_KMGW_inf_vs_un.tsv",sep = "\t",row.names = 1)

dea$Change <- "Other"
dea$Change[grep("^TLR",row.names(dea))] <- "TLRs"

dea$Change <- factor(dea$Change, levels = c("TLRs","Other"))

dea <- dea[order(-dea$padj),]
dea <- rbind(dea[grep("^TLR",row.names(dea),invert = T),], dea[grep("^TLR",row.names(dea)),] ) # TLR first

dea$txt <- NA
keep <- row.names(dea)[grepl("^TLR",row.names(dea)) & abs(dea$log2FoldChange) >1]
dea[keep,"txt"] <- keep

dea$size <- 2
dea[keep,"size"] <- 3
dea$shape <- 1
dea[keep,"shape"] <- 19

p1 <- ggplot(dea,aes(log2FoldChange,-log10(padj),color=Change,shape=shape))+scale_shape_identity()+
      geom_point(size=dea$size)+theme_bw()+
      #scale_color_manual(values=c("#FF000020","#0000FF20","#000000","#80808020"))+
      scale_color_manual(values=c("#000000","#a0a0a0"))+
      theme(panel.grid = element_blank(),axis.title = element_text(size=15),aspect.ratio = 1,
            plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5),
            axis.text = element_text(size = 13), legend.position = "none",
            legend.text = element_text(size = 13),legend.title = element_blank() )+
      geom_vline(xintercept = c(-1,1),linetype="dashed")+geom_hline(yintercept = -log10(0.05),linetype="dashed")+
      scale_x_continuous(limits = c(-20,20),breaks = 5*(-4:4),oob = scales::oob_squish)+ #axis.text = element_text(size = 13),
      scale_y_continuous(limits = c(0,50),oob = scales::oob_squish)+
      guides(color = guide_legend(override.aes = list(size = 4,alpha=1) ) )+
      ggtitle("Differential Expression Analysis","RNASeq, Infected vs Uninfected")+ 
      ggrepel::geom_text_repel(aes(label=dea$txt),nudge_x = -5.5,show.legend = FALSE)
print(p1)
#save
########### Boxplot ##########
tsv <- read.csv("C:/Users/Aris/Desktop/GSEA/Inputs/KMGW_B_inf_un.gct", sep="\t",skip = 2)[,-2]
toplot <- reshape::melt(tsv[grep("^TLR",tsv$NAME),])

colnames(toplot)[1:3]<-c("Gene","Sample","TPM")

#cp <- coord_polar(theta = "x"); cp$is_free <- function() TRUE

toplot$Sample <- factor(toplot$Sample,levels = unique(toplot$Sample))
toplot$Condition <- "Uninfected" ; toplot$Condition[grep("^T",toplot$Sample)] <- "Infected"
toplot$Condition <- factor(toplot$Condition,levels = c("Uninfected","Infected"))

toplot$Gene <- factor(toplot$Gene, levels = unique(toplot$Gene)[c(2,5,6,8,4,7,3,1)])


p2 <- ggplot(toplot, aes(Gene,TPM))+
        geom_boxplot(aes(fill=Condition),outlier.shape = NA)+
        scale_fill_manual(values=c("#0000FF40","#FF000040"))+
        geom_point(aes(color=Condition),show.legend = F,position = position_jitterdodge(jitter.width =0.15,dodge.width = 0.75))+
        scale_color_manual(values=c("#0000FF","#FF0000"))+
        xlab("Gene")+
        theme_bw()+
        theme(plot.title = element_text(hjust=0.5),plot.subtitle = element_text(hjust=0.5),panel.grid = element_blank(),
              aspect.ratio = 1,strip.text = element_text(size = 12),
              legend.text = element_text(size = 15),legend.title = element_blank(),
              axis.ticks.x=element_blank(),axis.title = element_text(size=15),axis.title.x = element_blank(),
              legend.position = "bottom",
              axis.text = element_text(size = 13,colour = "#000000"))+
        guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
        geom_vline(xintercept = 0.5+1:(length(unique(toplot$Gene)) -1 ))
print(p2)

library(ggpubr) ; library(cowplot)
ggarrange(plotlist = align_plots(p1, p2, axis = "tblr", align = "vh"))

## manual rsad2 plot uuuuh ##
lps <- openxlsx::read.xlsx("C:/Users/Aris/Downloads/Aris raw data LPS.xlsx",sheet=3)
#lps$value <- log10(lps$`2^-deltadeltaCt`)
lps$value <- lps$`2^-deltadeltaCt`
lps$cell <- c(rep("BL3",16),rep("TBL3",16))

lps$Sample <- gsub("BL3|TBL3|-|LPS","",lps$Sample)

lps$Sample <- factor(lps$Sample, levels =unique(lps$Sample))
lps$g <- rep(1:4,8)


lps <-na.omit(lps)

p1 <- ggplot(lps,aes(Sample,value))+
        geom_boxplot(outlier.alpha = 0)+
        geom_point(aes(group=g,color=cell),position = position_dodge(width = 0.25),size=3)+
        theme_bw()+
        scale_color_manual(values = c("#0000FF","#FF0000"))+
        theme(panel.grid = element_blank(),
              strip.text = element_text(size = 18,face="bold"),
              legend.text = element_text(size = 18),legend.title = element_blank(),
              axis.ticks.x=element_blank(),axis.title = element_text(size=15),
              legend.position = "none",strip.background = element_blank(),
              axis.text = element_text(size = 15,colour = "#000000")
              )+
        facet_wrap("cell",scales = "free_x")+
        ylab("RSAD2 mRNA expression\nFold Change")+
        xlab("LPS Concentrations (μg/mL)")+
        scale_y_continuous(breaks = 1:10)
print(p1)

####### Parva Analysis #######
tsv <- openxlsx::read.xlsx("C:/Users/Aris/Downloads/41598_2024_59197_MOESM4_ESM.xlsx")[,c(2,11:22)]
tsv <- tsv[!duplicated(tsv$Symbol),] ; row.names(tsv) <- tsv$Symbol ; #tsv <- tsv[,-1]

toplot <- reshape::melt(tsv[grep("^TLR",row.names(tsv)),])
colnames(toplot)[1:3]<-c("Gene","Sample","TPM")

toplot$Sample <- factor(toplot$Sample,levels = unique(toplot$Sample))
toplot$Condition <- "Uninfected" ; toplot$Condition[!grepl("Control",toplot$Sample)] <- "Infected"
toplot$Condition <- factor(toplot$Condition,levels = c("Uninfected","Infected"))

p2 <- ggplot(toplot, aes(Gene,TPM))+
  geom_boxplot(aes(fill=Condition),outlier.shape = NA)+
  scale_fill_manual(values=c("#0000FF40","#FF000040"))+
  geom_point(aes(color=Condition),show.legend = F,position = position_jitterdodge(jitter.width =0.15,dodge.width = 0.75))+
  scale_color_manual(values=c("#0000FF","#FF0000"))+
  xlab("Gene")+ylab("Normalized Counts")+
  theme_bw()+
  theme(plot.title = ggtext::element_markdown(hjust=0.5,face="bold",size = 12),plot.subtitle = element_text(hjust=0.5),panel.grid = element_blank(),
        aspect.ratio = 1,strip.text = element_text(size = 12),
        legend.text = element_text(size = 15),legend.title = element_blank(),
        axis.ticks.x=element_blank(),axis.title = element_text(size=15),axis.title.x = element_blank(),
        legend.position = "bottom",
        axis.text = element_text(size = 13,colour = "#000000"))+
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  geom_vline(xintercept = 0.5+1:(length(unique(toplot$Gene)) -1 ))+ggtitle("*Theileria parva* infection (7 days)")
print(p2)

library(cowplot)
plot_grid(plotlist = align_plots(p1,NULL,p2,align = "h",axis = "tb"),rel_widths = c(1,0,1),nrow = 1)

#### Figure 3 - GBPs ####
#########################
dea <- read.table("dea_TBL3_M_BL3_M.tsv")
dea$col <- "Other"
dea$col[grep("^GBP|LOC512486$|LOC507055$|LOC511531$|LOC513659$|LOC112445995$|LOC781710",row.names(dea),F)] <- "GBPs"

keep <- dea$col!="Other" & abs(dea$log2FoldChange) >1
dea$size <- 2
dea[keep,"size"] <- 3
dea$shape <- 1
dea[keep,"shape"] <- 19
dea$txt <- NA
dea$txt[keep] <- row.names(dea)[keep]
dea <- dea[order(dea$size),]
dea$padj[dea$padj==0] <- min(dea$padj[dea$padj!=0]) #cap

locs <- data.frame( LOCs = c("LOC512486","LOC781710","LOC513659","GBP5","LOC507055","GBP4","LOC112445995","LOC511531","GBP6"),
                    Name = c("GBP1","GBP2","GBP6_1","GBP5","GBP4_1","GBP4","GBP7-like","GBP1_1","GBP6"))
dea$lab <- locs[match(dea$txt,locs[,1]),2]




p1 <- ggplot(dea,aes(log2FoldChange,-log10(padj),color=col,shape=shape))+scale_shape_identity()+
  geom_point(size=dea$size)+theme_bw()+
  #scale_color_manual(values=c("#FF000020","#0000FF20","#000000","#80808020"))+
  scale_color_manual(values=c("#000000","#a0a0a0"))+
  theme(panel.grid = element_blank(),axis.title = element_text(size=15),aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5,size = 18),
        plot.subtitle = element_text(hjust = 0.5,size = 15),
        axis.text = element_text(size = 13),legend.position = "none",
        legend.text = element_text(size = 13),legend.title = element_blank() )+
  geom_vline(xintercept = c(-1,1),linetype="dashed")+geom_hline(yintercept = -log10(0.05),linetype="dashed")+
  scale_x_continuous(limits = c(-20,20),breaks = 5*(-4:4),oob = scales::oob_squish)+ #axis.text = element_text(size = 13),
  scale_y_continuous(limits = c(0,310),oob = scales::oob_squish)+
  guides(color = guide_legend(override.aes = list(size = 4) ) )+
  ggtitle("RNASeq, TBL3 vs BL3 (Villares et al.)")+ 
ggrepel::geom_text_repel(aes(label=dea$lab),nudge_x = -3.5,show.legend = FALSE)
print(p1)

#library(ggpubr) ; library(cowplot)
#ggarrange(plotlist = align_plots(p3, p2, axis = "tblr", align = "hv"),ncol = 1)

## protein GBPs
#dea2 <- read.csv("dea_prot_inf_un.tsv",sep="\t")
#dea2 <- dea2[!(dea2$Gene %in% readLines("TAgenes.txt")),]
dea2 <- readxl::read_xlsx("JD_proteomic_analysis.xlsx","Sheet1",col_types =  c("text","text","numeric","numeric","numeric","numeric"))
colnames(dea2)[c(1,5,6)] <- c("Gene","log2FoldChange","NegLog10pval")

dea2$col <- "Other"
dea2$col[grep("^GBP|LOC512486$|LOC507055$|LOC511531$|LOC513659$|LOC112445995$|LOC781710",dea2$Gene,F)] <- "GBPs"
dea2 <- dea2[!is.na(dea2$NegLog10pval),] #GBP4 is unique

keep <- dea2$col!="Other" 
dea2$size <- 2
dea2[keep,"size"] <- 3
dea2$shape <- 1
dea2[keep,"shape"] <- 19
dea2$txt <- NA
dea2$txt[keep] <- dea2$Gene[keep]
dea2 <- dea2[order(dea2$size),]
#dea2$padj[dea2$padj==0] <- min(dea2$padj[dea2$padj!=0]) #cap

locs <- data.frame( LOCs = c("LOC512486","GBP2","LOC781710","LOC513659","GBP5","LOC507055","GBP4","LOC112445995","LOC511531","GBP6"),
                    Name = c("GBP1","GBP2","GBP2","GBP6_1","GBP5","GBP4_1","GBP4","GBP7-like","GBP1_1","GBP6"))
dea2$lab <- locs[match(dea2$txt,locs[,1]),2]


p2 <- ggplot(dea2,aes(log2FoldChange,NegLog10pval,color=col,shape=shape))+scale_shape_identity()+
  geom_point(size=dea2$size)+theme_bw()+
  #scale_color_manual(values=c("#FF000020","#0000FF20","#000000","#80808020"))+
  scale_color_manual(values=c("#000000","#a0a0a0"))+ylab("-log10(pval)")+
  theme(panel.grid = element_blank(),axis.title = element_text(size=15),aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5,size = 18),
        plot.subtitle = element_text(hjust = 0.5,size = 15),
        axis.text = element_text(size = 13),legend.position = "none",
        legend.text = element_text(size = 13),legend.title = element_blank() )+
  geom_vline(xintercept = c(-1,1),linetype="dashed")+geom_hline(yintercept = -log10(0.05),linetype="dashed")+
  scale_x_continuous(limits = c(-10,10),breaks = 5*(-2:2),oob = scales::oob_squish)+ #axis.text = element_text(size = 13),
  scale_y_continuous(limits = c(0,7.5),oob = scales::oob_squish)+
  guides(color = guide_legend(override.aes = list(size = 4) ) )+
  ggtitle("Proteome, TBL3 vs BL3")+ 
  ggrepel::geom_text_repel(aes(label=dea2$lab),nudge_x = -3.5,show.legend = FALSE)+
  geom_rect(xmin = -12,xmax=-5,ymin=6.5,ymax=10,fill="white",lty="dashed")+
  annotate("text",-8,7.5,label="BL3 Unique",fontface="bold")+
  annotate("text",-8,6.87,label="GBP4")

print(p2)
library(cowplot)
plot_grid(plotlist = align_plots(p1,NULL,p2,axis = "h",align = "tb"),nrow = 1,rel_widths = c(1,0,1))

###################
tsv <- openxlsx::read.xlsx("C:/Users/Aris/Downloads/41598_2024_59197_MOESM4_ESM.xlsx")[,c(2,11:22)]
tsv <- tsv[!duplicated(tsv$Symbol),] ; row.names(tsv) <- tsv$Symbol ; #tsv <- tsv[,-1]
toplot <- reshape::melt(tsv[grep("^GBP|LOC512486$|LOC507055$|LOC511531$|LOC513659$|LOC112445995$|LOC781710",row.names(tsv)),])
colnames(toplot)[1:3]<-c("Gene","Sample","TPM")

toplot$Sample <- factor(toplot$Sample,levels = unique(toplot$Sample))
toplot$Condition <- "Uninfected" ; toplot$Condition[!grepl("Control",toplot$Sample)] <- "Infected"
toplot$Condition <- factor(toplot$Condition,levels = c("Uninfected","Infected"))


#toplot$Gene <- factor(toplot$Gene, levels = unique(toplot$Gene)[c(2,5,6,8,4,7,3,1)])

p3 <- ggplot(toplot, aes(Gene,TPM))+
  geom_boxplot(aes(fill=Condition),outlier.shape = NA)+
  scale_fill_manual(values=c("#0000FF40","#FF000040"))+
  geom_point(aes(color=Condition),show.legend = F,position = position_jitterdodge(jitter.width =0.15,dodge.width = 0.75))+
  scale_color_manual(values=c("#0000FF","#FF0000"))+
  xlab("Gene")+ylab("Normalized Counts")+scale_y_continuous(limits = c(0,300),oob = scales::oob_squish)+
  theme_bw()+
  theme(plot.title = element_text(hjust=0.5),plot.subtitle = element_text(hjust=0.5),panel.grid = element_blank(),
        aspect.ratio = 1,strip.text = element_text(size = 12),
        legend.text = element_text(size = 15),legend.title = element_blank(),
        axis.ticks.x=element_blank(),axis.title = element_text(size=15),axis.title.x = element_blank(),
        legend.position = "bottom",
        axis.text = element_text(size = 13,colour = "#000000"))+
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  geom_vline(xintercept = 0.5+1:(length(unique(toplot$Gene)) -1 ))
print(p3)


########## FIGURE 4 ##########
##########  Heatmap ##########
tsv <- read.csv("C:/Users/Aris/Desktop/GSEA/Inputs/KMGW_B_inf_un.gct", sep="\t",skip = 2)[,-2]
toplot <- reshape::melt(tsv[grep("^BOLA",tsv$NAME,ignore.case = T),])

colnames(toplot)[1:3]<-c("Gene","Sample","TPM")

toplot$Sample <- factor(toplot$Sample,levels = rev(unique(toplot$Sample)))
toplot$Gene <- factor(toplot$Gene,levels = unique(c(grep("R",toplot$Gene,value = T),grep("R",toplot$Gene,value = T,invert = T) )))
toplot$Condition <- "Uninfected" ; toplot$Condition[grep("^T",toplot$Sample)] <- "Infected"
toplot$Condition <- factor(toplot$Condition,levels = c("Uninfected","Infected"))

toplot$MHC<- "Class I"
toplot$MHC[grep("-D",toplot$Gene)]<- "Class II"

cp <- coord_polar(theta = "x"); cp$is_free <- function() TRUE


p1 <- ggplot(toplot, aes(Sample,reorder(Gene,TPM),fill=log10(TPM+1)))+
        geom_tile()+
        viridis::scale_fill_viridis(option = "inferno",direction = -1)+
        theme_bw()+ylab("Gene")+xlab("")+
        theme(plot.title = element_text(hjust=0.5),plot.subtitle = element_text(hjust=0.5),panel.grid = element_blank(),
              aspect.ratio = 1,strip.background = element_blank(),
              strip.text = element_text(size = 12),strip.placement = "outside",
              legend.text = element_text(size = 13),
              panel.border = element_blank(),
              axis.ticks.x=element_blank(),axis.title = element_text(size=15),
              axis.ticks.length.y=unit(-3.3, units = "cm"),axis.title.y = element_blank(),
              axis.ticks.y = element_line(linetype = "dotted",linewidth = 0.35,colour="#ffffff"),
              axis.text.x = element_text(face = "bold",size = 13,colour = "#000000"), #element_blank(),#element_text(angle = 90,hjust = 1),
              axis.text.y = element_text(hjust = 0),
              plot.background = element_blank(),
              axis.text = element_text(colour = "#000000"))+
        scale_x_discrete(labels=c(rep("",12),"\nInfected    ",rep("",4),"\nUninfected ",rep("",12)))+
        facet_wrap("MHC",nrow = 2,scales = "free_y",strip.position = "left")+
        cp
print(p1)

dea <- read.table("dea_TBL3_M_BL3_M.tsv")
dea$Gene <- row.names(dea)
dea$col <- "Other"
dea$col[dea$Gene == "MMP9"] <- "MMP9"
dea$col[grep("BOLA",dea$Gene,F)] <- "MHC Class I"
dea$col[grep("BOLA-D|LOC100848815$",dea$Gene)] <- "MHC Class II"
dea$size <- 1
dea$size[dea$col!="Other" & abs(dea$log2FoldChange) >1] <- 4
dea$shape <- 1
dea$shape[dea$col!="Other" & abs(dea$log2FoldChange) >1] <- 19
dea$txt <- NA
dea$txt[dea$col!="Other"] <- dea$Gene[dea$col!="Other"]

dea$padj[dea$padj==0] <- min(dea$padj[dea$padj!=0])
dea <- dea[order(dea$size),]

p2 <- ggplot(dea,aes(log2FoldChange,-log10(padj),color=col,shape=shape))+scale_shape_identity()+
  geom_point(size=dea$size)+theme_bw()+
  scale_color_manual(values=c("#207F00","#A000F0","#000000","#a0a0a0"))+
  #scale_color_manual(values=c("#000000","#80808020"))+
  theme(panel.grid = element_blank(),axis.title = element_text(size=15),aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5,size = 18),
        plot.subtitle = element_text(hjust = 0.5,size = 15),
        axis.text = element_text(size = 13),#legend.position = "none",
        legend.text = element_text(size = 13),legend.title = element_blank() )+
  geom_vline(xintercept = c(-1,1),linetype="dashed")+geom_hline(yintercept = -log10(0.05),linetype="dashed")+
  scale_x_continuous(limits = c(-20,20),breaks = 5*(-4:4),oob = scales::oob_squish)+ #axis.text = element_text(size = 13),
  scale_y_continuous(limits = c(0,310),oob = scales::oob_squish)+
  guides(color = guide_legend(override.aes = list(size = 4) ) )+
  ggtitle("RNASeq, TBL3 vs BL3 (Villares et al.)")#+ 
#ggrepel::geom_text_repel(aes(label=dea$txt),nudge_x = 3.5,show.legend = FALSE)
print(p2)


######  Proteomics DEA  ######
#library(DESeq2)
#psv <- read.csv("C:/Users/Aris/Documents/Immune Evasion Paper/Theileria_TBL3-BL3_ProtDataFIXED.tsv", sep="\t",row.names = 1)
#todea <- data.frame(sapply(psv[,grep("^Abundance_",colnames(psv))[1:6]   ],as.numeric),row.names = row.names(psv))
#todea <- round(todea/min(na.omit(unlist(todea))))

#todea[is.na(todea)] <- min(na.omit(unlist(todea)))



#colD<- data.frame(
#  sample=colnames(todea),
#  condition=c(rep("BL3",3),rep("TBL3",3))
#)

## create the DESeqDataSet
#dds <- DESeqDataSetFromMatrix(todea, colData = colD, design = ~condition)

# generate normalized counts
#dds <- estimateSizeFactors(dds)

## perform DEA
#dds <- DESeq(dds)


###dea infected vs uninfected
#dea2 <- na.omit(cbind(as.data.frame(results(dds, contrast=c("condition","TBL3","BL3"))),psv$GeneId)) ; colnames(dea2)[7]<-"Gene"
#write.table(dea.sorted <- dea[order(dea$padj, -abs(dea$log2FoldChange), decreasing = FALSE), ],
#            "Immune Evasion Paper/dea_prot_inf_un.tsv",col.names = NA,quote = F,sep = "\t")  # sort the table: ascending of padj then descending of absolute valued of logFC

dea2 <- readxl::read_xlsx("C:/Users/Aris/Documents/Immune Evasion Paper/JD_proteomic_analysis.xlsx","Sheet1",col_types =  c("text","text","numeric","numeric","numeric","numeric"))
colnames(dea2)[c(1,5,6)] <- c("Gene","log2FoldChange","NegLog10pval")
dea2$Gene <- toupper(dea2$Gene)

dea2$col <- "Other"
dea2$col[dea2$Gene == "MMP9"] <- "MMP9"
dea2$col[grep("BOLA",dea2$Gene,F)] <- "MHC Class I"
dea2$col[grep("BOLA-D|LOC100848815$",dea2$Gene)] <- "MHC Class II"
dea2$size <- 1
dea2$size[dea2$col!="Other" ] <- 4
dea2$shape <- 1
dea2$shape[dea2$col!="Other" ] <- 19
dea2$txt <- NA
dea2$txt[dea2$col!="Other"] <- dea2$Gene[dea2$col!="Other"]

#dea2$NegLog10pval[dea2$NegLog10pval==0] <- min(dea2$NegLog10pval[dea2$NegLog10pval!=0])
dea2 <- dea2[order(dea2$size),]

p3 <- ggplot(dea2,aes(log2FoldChange,NegLog10pval,color=col,shape=shape))+scale_shape_identity()+
        geom_point(size=dea2$size)+theme_bw()+ylab("-log10(pval)")+
        scale_color_manual(values=c("#207F00","#A000F0","#000000","#a0a0a0"))+
        #scale_color_manual(values=c("#000000","#80808020"))+
        theme(panel.grid = element_blank(),axis.title = element_text(size=15),aspect.ratio = 1,
              plot.title = element_text(hjust = 0.5,size = 18),
              plot.subtitle = element_text(hjust = 0.5,size = 15),
              axis.text = element_text(size = 13),#legend.position = "none",
              legend.text = element_text(size = 13),legend.title = element_blank() )+
        geom_vline(xintercept = c(-1,1),linetype="dashed")+geom_hline(yintercept = -log10(0.05),linetype="dashed")+
        scale_x_continuous(limits = c(-10,10),breaks = 5*(-2:2),oob = scales::oob_squish)+ #axis.text = element_text(size = 13),
        scale_y_continuous(limits = c(0,7.5),oob = scales::oob_squish)+
        guides(color = guide_legend(override.aes = list(size = 4) ) )+
        ggrepel::geom_text_repel(aes(label=dea2$txt),show.legend = FALSE,nudge_x = 4.5*sign(dea2$log2FoldChange),force = 0.5)+
        geom_rect(xmin = -12,xmax=-5,ymin=5.5,ymax=10,fill="white",color="black",lty="dashed",show.legend = F)+
        annotate("text",-8,7.5,label="BL3 Unique",fontface="bold")+
        annotate("text",-8,6.5,label="BOLA-DQB\nBOLA-DRB2\nBOLA-DRB3",color = "#A000F0")+
        ggtitle("Proteome, Infected vs Uninfected")#+ 
        #ggrepel::geom_text_repel(aes(label=dea2$txt),nudge_x = 3.5,show.legend = FALSE)
print(p3)



#"RNASeq, TBL3 vs BL3 (Villares et al.)"
library(ggpubr) ; library(cowplot)
ggarrange(plotlist = align_plots(p2, p3, axis = "tblr", align = "hv"),ncol = 1)



########## FIGURE 6 ##########
#########  CIITA TPM #########
tsv1 <- read.table("TPM.tsv")

toplot <- reshape::melt(cbind(tsv1,row.names(tsv1))[grep("^CIITA|^CD74$",row.names(tsv1)),])

colnames(toplot) <- c("Gene","Sample","TPM")
toplot$Condition <- gsub("_\\d$","",toplot$Sample)

toplot$Dataset <- "Cheeseman et al."
toplot$Dataset[grep("M|2646",toplot$Sample)] <- "Villares et al."

toplot$Condition <- gsub("_K$|_M$","",toplot$Condition)
toplot$Condition <- factor(toplot$Condition, levels = c("BL3", "BL3_Bup","BL3_2646", "TBL3", "TBL3_Bup", "TBL3_2646"))

#toplot$Gene <- factor(toplot$Gene, levels = c("CIITA","RFX5","NFYB","CREB1","CD74"))
toplot$Gene <- factor(toplot$Gene, levels = c("CIITA","CD74"))


ggplot(toplot, aes(Condition,log10(TPM+1)))+xlab("")+ylab("Gene Expression (Log10)")+
  geom_boxplot(aes(fill=Condition),outlier.shape = NA)+
  scale_fill_manual(values=c("#0000FF40","#0000FF40","#0000FF40","#FF000040","#FF000040","#FF000040"))+
geom_point(aes(color=Condition),position = position_jitterdodge(jitter.width =0.15,dodge.width = 0.75))+
  scale_color_manual(values=c("#0000FF","#0000FF40","#0000FF20","#FF0000","#FF000040","#FF000020"))+
  theme_bw()+
  theme(plot.title = element_text(hjust=0.5),plot.subtitle = element_text(hjust=0.5),panel.grid = element_blank(),
        aspect.ratio = 1,
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 13),legend.title = element_blank(),legend.position = "none",
        axis.ticks.x=element_blank(),axis.title = element_text(size=15),
        axis.text.x = element_text(size = 13,colour = "#000000",angle = 90,hjust=1),
        axis.text.y = element_text(colour = "#000000"))+
  guides(color = guide_legend(override.aes = list(size = 5)))+
  facet_wrap("Dataset~Gene",scales = "free_x",nrow = 2)

########## FIGURE 7 ##########
#####  Reversibility TPM #####
fp <- list.files(".","counts.tsv",full.names = T)
names(fp) <- list.files(".","counts.tsv")
tagenes <- readLines("../TAgenes.txt")

tsv <- Reduce(cbind,lapply(names(fp),function(x){
  l <- read.csv(fp[x],sep = "\t",row.names = 1)
  colnames(l) <- paste(gsub("_counts\\.tsv","",x),"_",1:3,sep="")
  l 
}))


#### TPM ####
#time to normalize for genesize
library(GenomicFeatures)
txdb <- loadDb("./Bt_Ta_fused.txdb") 
tx_by_gene <- transcriptsBy(txdb, by="gene")
tsize <- max(abs(end(tx_by_gene)-start(tx_by_gene))) #keep largest transcript per gene

tsize<-tsize[row.names(tsv)]

tsv1<-tsv[names(tsize),] #put counts in same order as in transcript sizes
tsv1 <- sweep(tsv1, 1, tsize, "/") #normalize genesize
tsv1<-na.omit(tsv1)
tsv1<- sweep(tsv1, 2, colSums(tsv1), "/") #normalize for total transcripts
colSums(tsv1) # check new sizes if same

tsv1<-tsv1*1e6 #multiply by million to get TPM
#write.table(tsv1,".TPM.tsv",F,F,"\t",col.names = NA)
tsv1 <- tsv1[rowSums(tsv1) != 0,]

toplot <- reshape::melt(cbind(tsv1,row.names(tsv1))[grep("^TLR|^BOLA-D|^GBP|LOC512486$|LOC507055$|LOC511531$|LOC513659$|LOC112445995$|LOC781710",row.names(tsv1)),])
#toplot <- reshape::melt(cbind(tsv1,row.names(tsv1))[grep("^NLR|^NOD",row.names(tsv1)),])

colnames(toplot) <- c("Gene","Sample","TPM")
toplot$Condition <- gsub("_\\d$","",toplot$Sample)

toplot$Condition <- factor(toplot$Condition, levels = c("BL3_K", "BL3_Bup", "TBL3_K", "TBL3_Bup", "BL3_M","BL3_2646", "TBL3_M", "TBL3_2646"))
toplot$Group <- "Innate immunity"
toplot$Group[grep("BOLA-D",toplot$Gene)] <- "Acquired immunity"
toplot$Group[grep("^GBP|LOC512486$|LOC507055$|LOC511531$|LOC513659$|LOC112445995$|LOC781710",toplot$Gene)] <- "Inflammasome activation"
toplot$Group <- factor(toplot$Group,levels = c("Innate immunity","Inflammasome activation","Acquired immunity"))
toplot$Dataset <- "Cheeseman et al."
toplot$Dataset[grep("M|2646",toplot$Sample)] <- "Villares et al."

labs <- rep("",24 ); labs[c(2 + 3*0:7)] <- c("BL3", "BL3_Bup", "TBL3", "TBL3_Bup", "BL3","BL3_2646", "TBL3", "TBL3_2646")

locs <- data.frame( LOCs = c("LOC512486","LOC781710","LOC513659","GBP5","LOC507055","GBP4","LOC112445995","LOC511531","GBP6",grep("^TLR|^BOLA-D",toplot$Gene,value = T)),
                    Name = c("GBP1","GBP2","GBP6_1","GBP5","GBP4_1","GBP4","GBP7-like","GBP1_1","GBP6",grep("^TLR|^BOLA-D",toplot$Gene,value = T)))
toplot$Gene <- locs[match(toplot$Gene,locs[,1]),2]


ggplot(toplot, aes(reorder(Sample,as.numeric(Condition)),reorder(Gene,TPM),fill=log10(TPM+1)))+
  geom_tile(show.legend = TRUE)+
  theme_bw()+ylab("Gene")+
  viridis::scale_fill_viridis(option = "inferno",direction = -1)+ #inferno/mako
  ggtitle("Reversibility of Immune Genes","RNASeq, Treated vs Untreated") + 
  theme(plot.title = element_text(hjust=0.5),plot.subtitle = element_text(hjust=0.5),panel.grid = element_blank(),
        aspect.ratio = 1,strip.text = element_text(size = 12),
        axis.title = element_blank(),
        axis.ticks.x=element_blank(),strip.background = element_blank(),
        axis.text.x = element_text(size = 10, colour = "#000000",face="bold"),
        axis.text.y = element_text(size = 10,colour = "#000000",hjust=0))+
  scale_x_discrete(breaks=levels(toplot$Sample)[c(7:9,4:6,19:21,16:18,10:12,1:3,22:24,13:15)],labels=labs,expand = expansion(mult = -0.1))+ # reorganize labels coz R hates me
  scale_y_discrete(expand = expansion(mult = -0.1))+
  facet_wrap("Group~Dataset",scales = "free",ncol=2)+
  geom_vline(xintercept = c(0.5 + 3*1:3),lty="dashed")

#######  Supplementary #######
nms <- locs[!duplicated(locs$Name),] ; nms <- nms[rev(order(nms$LOCs)),]
nms <- rbind(nms,data.frame("LOCs"=c("CIITA","CD74"), "Name"=c("CIITA","CD74")))
gns <- genes(txdb)
gns <- gns[gns@ranges@NAMES %in% nms$LOCs]

library(rtracklayer)
export.bed(gns,con='granges.bed') 

