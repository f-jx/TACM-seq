.libPaths("/data/R02/nieyg/ori/biosoft/conda/envs/r4-base/lib/R/library")

library(BSgenome)
library(Matrix)
library(Signac)
library(Seurat)
library(BSgenome.Mmusculus.UCSC.mm10)
library(tidyverse)
library(cowplot)
library(patchwork)
library(sctransform)
library(scater)
library(celda)
library(hdf5r)

OSN_mm<-readRDS("./02_All_celltype/WNN_OSN_integrated_all_celltype.rds")

Neutrophils<- subset(OSN_mm,idents="Neutrophils")

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(Neutrophils) <- "ATAC"
########## 这里没有重新callpeak，直接用原来的矩阵数据重新UMAP并分群。
Neutrophils <- RunUMAP(Neutrophils, reduction = 'integrated_LSI', dims = 2:40, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

# tune the resolution
Neutrophils <- FindNeighbors(object = Neutrophils, reduction = 'integrated_LSI', dims = 2:40)

Neutrophils <- FindClusters(Neutrophils, resolution =1, algorithm = 3, verbose = FALSE)

my47colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
         '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
         '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
         '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', '#116530', '#678983',
         '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76')
pdf("./03_Neutrophils/Neutrophils_cluster_ATAC.pdf",width=15,height=5)
###cluster
p1 <- DimPlot(Neutrophils, cols=my47colors, reduction = "umap.atac", group.by = "seurat_clusters", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("seurat_clusters")
p2 <- DimPlot(Neutrophils, cols=my47colors, reduction = "umap.atac", group.by = "orig.ident", shuffle=TRUE,label = F, label.size = 5, repel = TRUE) + ggtitle("sample")
p3 <- FeaturePlot(Neutrophils, reduction = 'umap.atac',features = "atac_S_M.pulmonis_ReadsPercentage", ncol = 1)
p1 + p2 +p3 
dev.off()


########## 这个思路不错，可以调整resolution后，将小群按主观意愿合并并重命名。
Idents(Neutrophils)<-Neutrophils$seurat_clusters
#####further annotation########
Neutrophils <- RenameIdents(
  object = Neutrophils,
  '0' = 'C1',
  '1' = 'C2',
  '2' = 'C2',
  '3' = 'C3',
  '4' = 'C1',
  '5' = 'C3',
  '6' = 'C3',
  '7' = 'C4',
  '8' = 'C4'
  )
Neutrophils@meta.data$subtype<-Idents(Neutrophils)
table(Neutrophils$subtype,Neutrophils$orig.ident)

Neutrophils$subtype<-factor(Neutrophils$subtype,
  levels=c("C1","C2","C3","C4"
    ))
Idents(Neutrophils)<-Neutrophils$subtype;

# Dimplot with annotation
pdf("./03_Neutrophils/Neutrophils_cluster_ATAC_subtype.pdf",width=10,height=8)
###cluster
p1 <- DimPlot(Neutrophils, cols=my47colors[-1:-2], reduction = "umap.atac", group.by = "seurat_clusters", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("seurat_clusters")
p2 <- DimPlot(Neutrophils, cols=my47colors[-1:-2], reduction = "umap.atac", group.by = "subtype", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("subtype")
pp1<- (p1 + p2 )
p3 <- DimPlot(Neutrophils, cols=my47colors, reduction = "umap.atac", group.by = "orig.ident", shuffle=TRUE,label = F, label.size = 5, repel = TRUE) + ggtitle("sample")
p4 <- FeaturePlot(Neutrophils, reduction = 'umap.atac',features = "atac_S_M.pulmonis_ReadsPercentage", ncol = 1)
pp2<- (p3 + p4 )
pp1/pp2
dev.off()


# cell proportion distribution 
pdf("./03_Neutrophils/Neutrophils_subtype_proportion.pdf",width=3,height=5)
df <- as.data.frame(Neutrophils@meta.data)
df_ct <- df %>%
    group_by(orig.ident, subtype) %>%
    summarise(counts = n()) %>%
    mutate(cell_proportion = counts / sum(counts))

p <- ggplot(df_ct, aes(orig.ident, cell_proportion, fill=subtype)) + 
    geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = my47colors[-1:-2]) +
    theme_cowplot() +
    xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
p
dev.off()


pdf("./03_Neutrophils/Neutrophils_subtype_M.pulmonis_violinplot.pdf",width=4,height=8)

p1 <- VlnPlot(subset(Neutrophils,atac_S_M.pulmonis_ReadsPercentage>0), "atac_S_M.pulmonis_ReadsPercentage",add.noise=FALSE, pt.size = 0,col=my47colors[-1:-2],group.by = "subtype") +
geom_boxplot(width=0.1,outlier.size=0.01, position = position_dodge(0.9))

p2 <- VlnPlot(Neutrophils, "atac_S_M.pulmonis_ReadsPercentage", pt.size = 0.01,col=my47colors[-1:-2],group.by = "subtype") +
geom_boxplot(width=0.1,outlier.size=0.01, position = position_dodge(0.9))
p1/p2
dev.off()

########## 做差异分析前会重新Callpeak，并按新的矩阵计算DEP
# Find DEP for C1,2,3,4
# recall peak 
library(BSgenome.Mmusculus.UCSC.mm10)
# ORN recall peak for subcluster 
DefaultAssay(Neutrophils)<-"ATAC"
peak<-CallPeaks(
       Neutrophils,
       group.by = "subtype",
       macs2.path = "/public/home/nieyg/biosoft/conda/bin/macs2",
       broad = FALSE,
       format = "BED",
       fragment.tempdir = tempdir(),
       effective.genome.size = 1.87e+09, #mm10 size
       outdir="./03_Neutrophils/",
       combine.peaks=TRUE
)
macs2_counts <- FeatureMatrix(
     fragments = Fragments(Neutrophils),
     features = peak,
     cells = colnames(Neutrophils)
     )     
Neutrophils[["peaks_All_cluster"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  #fragments = fragpath[[i]],
  fragments = Fragments(Neutrophils),
  annotation = Annotation(Neutrophils)
)










