library(Seurat)
HC_drg1_0.data <- Read10X("../rawdata/HC-drg1/0/", gene.column=1)
HC_drg1_0 <- CreateSeuratObject(counts = HC_drg1_0.data, project = "HC_drg")

hist(colSums(HC_drg1_0),
     breaks = 150, main = "UMI count per cell",
     xlab = "UMI count per cell")
HC_drg1_0[["percent.mt"]] <- PercentageFeatureSet(HC_drg1_0, pattern = "^mt-")

VlnPlot(HC_drg1_0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.2)

FeatureScatter(HC_drg1_0, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(HC_drg1_0, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
HC_drg1_0 <- subset(HC_drg1_0, subset = nFeature_RNA > 100  & percent.mt < 5)
VlnPlot(HC_drg1_0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.2)

FeatureScatter(HC_drg1_0, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(HC_drg1_0, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

HC_drg1_0 <- NormalizeData(HC_drg1_0, normalization.method = "LogNormalize", scale.factor = 10000)

hist(colSums(HC_drg1_0$RNA@data),
     breaks = 150,
     main = "UMI count per cell after normalization",
     xlab = "UMI count per cell")

HC_drg1_0 <- FindVariableFeatures(HC_drg1_0, selection.method = "vst", nfeatures = 2000)
VariableFeaturePlot(HC_drg1_0)
HC_drg1_0 <- ScaleData(HC_drg1_0, features = rownames(HC_drg1_0))

HC_drg1_0 <- RunPCA(HC_drg1_0,verbose=FALSE)
ElbowPlot(HC_drg1_0, ndims= 40)
DimHeatmap(HC_drg1_0, dims = 1:10, cells = 100, balanced = TRUE)

HC_drg1_0 <- FindNeighbors(HC_drg1_0, dims = 1:10)
HC_drg1_0 <- FindClusters(HC_drg1_0, resolution = 0.5)

HC_drg1_0 <- RunUMAP(HC_drg1_0, dims = 1:10, n.neighbors = 20L, min.dist = 0.3,verbose=FALSE)
DimPlot(HC_drg1_0, reduction = "umap")

library(DoubletFinder)
sweep.res.list <- paramSweep_v3(HC_drg1_0, PCs = 1:10, sct = FALSE) 
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
p<-as.numeric(as.vector(bcmvn[bcmvn$MeanBC==max(bcmvn$MeanBC),]$pK))

homotypic.prop <- modelHomotypic(HC_drg1_0@meta.data$seurat_clusters)
nExp_poi <- round(Doubletrate*ncol(HC_drg1_0))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

HC_drg1_0 <- doubletFinder_v3(HC_drg1_0, PCs = 1:10, pN = 0.25, pK = p, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
colnames(HC_drg1_0@meta.data)[ncol(HC_drg1_0@meta.data)] = "doublet_info"

DimPlot(HC_drg1_0, reduction = "umap", group.by = "doublet_info")
HC_drg1_0 <- subset(HC_drg1_0,subset=doublet_info=="Singlet")
dim(HC_drg1_0)
save(HC_drg1_0, file = "HC_drgfiltered.rdata")

DRG_merged <- merge(HC_drg1_0, y = c(HC_drg1_1, con_drg2_0, con_drg2_1), 
                      add.cell.ids = c("HC_drg1_0", "HC_drg1_1", "con_drg2_0", "con_drg2_1"))
View(DRG_merged@meta.data)

DRG_merged <- NormalizeData(DRG_merged, normalization.method = "LogNormalize", scale.factor = 10000)
DRG_merged <- FindVariableFeatures(DRG_merged, selection.method = "vst", nfeatures = 2000,verbose=FALSE)
DRG_merged <- ScaleData(DRG_merged, features = rownames(DRG_merged))
DRG_merged <- RunPCA(DRG_merged,verbose=FALSE)
ElbowPlot(DRG_merged,ndims= 30)

DRG_merged <- FindNeighbors(DRG_merged, dims = 1:12, reduction = "pca")
DRG_merged <- FindClusters(DRG_merged, resolution = 0.5)
DRG_merged <- RunUMAP(DRG_merged, dims = 1:12, verbose=FALSE)
DRG_merged <- RunTSNE(DRG_merged, dims = 1:12, check_duplicates = FALSE)

library(clustree)
DRG_merged <- FindClusters(
  object = DRG_merged,
  resolution = c(seq(0,1.6,.2))
)
clustree(DRG_merged@meta.data, prefix = "RNA_snn_res.")

library(dplyr)
library(ggrepel)
library(ggthemes)
library(ggplot2)

DimPlot(DRG_merged, reduction = "umap", label=TRUE)+NoLegend()
DimPlot(DRG_merged, reduction = "tsne", label=TRUE)+NoLegend()

DRG_merged$orig.ident <- ifelse(grepl("^HC", DRG_merged$orig.ident), "HC_DRG", "con_DRG")
DimPlot(DRG_merged, reduction = "umap",group.by = "orig.ident", cols = c('con_DRG'='#6D91CB','HC_DRG'='#EFAFB4')) +
  theme(aspect.ratio = 1)

DefaultAssay(DRG_merged) <- "RNA"
DRG_merged_RNA_markers <- FindAllMarkers(DRG_merged, only.pos = TRUE, min.pct = 0.25,logfc.threshold = 0.1)
DRG_merged_all_markers <- FindAllMarkers(DRG_merged)
write.csv(DRG_merged_RNA_markers, "DRG_merged_markers.csv")
write.csv(DRG_merged_all_markers, "DRG_merged_all_markers.csv")

clu_marker <- read.csv("cluster_marker_per3.csv")
pdf("myplot.pdf", width = 4, height = 3)
DoHeatmap(DRG_merged,
          features = as.character(unique(clu_marker$gene)),
          group.by = "seurat_clusters",
          assay = 'RNA',
          angle = 0,
          draw.lines = F)+
  scale_fill_gradientn(colors = c("yellow","white","red"))
dev.off()

table(Idents(DRG_merged), DRG_merged$seurat_clusters)
seldrg_HC = DRG_merged[,DRG_merged@meta.data$orig.ident %in% "HC_DRG"]
table(Idents(seldrg_HC), seldrg_HC$seurat_clusters)
seldrg_con = DRG_merged[,DRG_merged@meta.data$orig.ident %in% "con_DRG"]
table(Idents(seldrg_con), seldrg_con$seurat_clusters)

load('DRG_merged.Rdata')
sample_table <- as.data.frame(table(DRG_merged@meta.data$orig.ident,DRG_merged@meta.data$seurat_clusters))
names(sample_table) <- c("con_HC","Cluster","CellNumber")
library(tidyverse)
library(patchwork)
colour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",  
         "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
         "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
         "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
plot_sample<-ggplot(sample_table,aes(x = con_HC, weight = CellNumber, fill = Cluster))+
  geom_bar(position="fill")+
  scale_fill_manual(values=colour) + 
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        axis.line.x = element_line(colour = "black") ,
        axis.line.y = element_line(colour = "black") ,
        plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size =16)
  )+labs(y="Percentage")+RotatedAxis()
pdf(file = 'DRG_per.pdf',width = 10, height = 10)
dev.off()

VlnPlot(DRG_merged, features = "TRPA1")                      
FeaturePlot(DRG_merged, features = "TRPA1", pt.size = 0.5, label = T)

seldrg_clu4 = DRG_merged[,DRG_merged@meta.data$seurat_clusters %in% c(4)]
seldrg_clu4_a1 <- subset(x = seldrg_clu4, subset = TRPA1 > 0)
marker_sel4_a1 <- FindMarkers(seldrg_clu4_a1, ident.1 = "con_DRG", ident.2 = "HC_DRG",group.by = 'orig.ident')
write.csv(marker_sel4_a1,file='con-hc-clu4a1.csv')

DEG_limma4 <- read.csv("con-hc-clu4a1.csv", header=T)
cut_off_pvalue = 0.05
cut_off_logFC = 1

DEG_limma4 <- 
  DEG_limma4 %>% 
  mutate(change = as.factor(ifelse(DEG_limma4$p_val < 0.05 & abs(DEG_limma4$avg_log2FC) > 1,
                                   ifelse(DEG_limma4$avg_log2FC > 1 ,'Up','Down'),'Stable')))
p <- ggplot(
  DEG_limma4, 
  aes(x = avg_log2FC, 
      y = -log10(p_val), 
      colour=change)) +
  geom_point(alpha=1, size=3.5) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
  
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_pvalue),lty=4,col="black",lwd=0.8) +
  
  labs(x="avg_log2FC",
       y="-log10 (p-value)")+
  theme_base()+
  
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank()
  )
gene <- read.csv(file = 'gene.csv')
gene
gene$geneList <- gene$genename
gene
data2 <- DEG_limma4 %>% left_join(gene, by = c("gene" = "genename"))
p + geom_text_repel(data = data2, aes(x = avg_log2FC, 
                                      y = -log10(p_val), 
                                      label = geneList),
                    size = 5,box.padding = unit(1, "lines"),
                    point.padding = unit(0.8, "lines"), 
                    segment.color = "black", 
                    max.overlaps =50,
                    show.legend = FALSE,
                    color="black"
)

BiocManager::install("org.Gg.eg.db")
library(org.Gg.eg.db)
library(clusterProfiler)
library(DOSE)
femarker <- read.csv("clu4_chayi.csv")
gene.df <- bitr(femarker$gene,fromType="SYMBOL",
	toType=c("ENSEMBL", "ENTREZID"), 
	OrgDb = org.Gg.eg.db)

gene <- gene.df$ENTREZID
ego_BP <- enrichGO(gene = gene,
                   OrgDb=org.Gg.eg.db,
                   keyType = "ENTREZID",
                   ont = "BP",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 1,
                   qvalueCutoff = 1,
                   readable = TRUE)
ego_result_BP <- as.data.frame(ego_BP)
write.csv(ego_result_BP,file = 'clu4_chayi_BP.csv')
pdf("3clu4_chayi_BP.pdf", width = 8, height = 10)
dotplot(ego_BP,title="EnrichmentGO_BP_markerdot")
dev.off()
