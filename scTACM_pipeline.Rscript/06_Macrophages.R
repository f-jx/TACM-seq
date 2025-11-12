
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

Macrophages<- subset(OSN_mm,idents="Macrophages")

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(Macrophages) <- "ATAC"
Macrophages <- RunUMAP(Macrophages, reduction = 'integrated_LSI', dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
Macrophages <- FindNeighbors(object = Macrophages, reduction = 'integrated_LSI', dims = 2:30)
Macrophages <- FindClusters(Macrophages, resolution =0.5, algorithm = 3, verbose = FALSE)

my47colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
         '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
         '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
         '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', '#116530', '#678983',
         '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76')
pdf("./04_Macrophages/Macrophages_cluster_ATAC.pdf",width=15,height=5)
###cluster
p1 <- DimPlot(Macrophages, cols=my47colors, reduction = "umap.atac", group.by = "seurat_clusters", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("seurat_clusters")
p2 <- DimPlot(Macrophages, cols=my47colors, reduction = "umap.atac", group.by = "orig.ident", shuffle=TRUE,label = F, label.size = 5, repel = TRUE) + ggtitle("sample")
p3 <- FeaturePlot(Macrophages, reduction = 'umap.atac',features = "atac_S_M.pulmonis_ReadsPercentage", ncol = 1)
p1 + p2 +p3 
dev.off()


Idents(Macrophages)<-Macrophages$seurat_clusters

Macrophages@meta.data$subtype<-Macrophages$seurat_clusters
table(Macrophages$subtype,Macrophages$orig.ident)
Idents(Macrophages)<-Macrophages$subtype;

# Dimplot with annotation
pdf("./04_Macrophages/Macrophages_cluster_ATAC_subtype.pdf",width=10,height=8)
###cluster
p1 <- DimPlot(Macrophages, cols=my47colors[-1:-6], reduction = "umap.atac", group.by = "seurat_clusters", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("seurat_clusters")
p2 <- DimPlot(Macrophages, cols=my47colors[-1:-6], reduction = "umap.atac", group.by = "subtype", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("subtype")
pp1<- (p1 + p2 )
p3 <- DimPlot(Macrophages, cols=my47colors, reduction = "umap.atac", group.by = "orig.ident", shuffle=TRUE,label = F, label.size = 5, repel = TRUE) + ggtitle("sample")
p4 <- FeaturePlot(Macrophages, reduction = 'umap.atac',features = "atac_S_M.pulmonis_ReadsPercentage", ncol = 1)
pp2<- (p3 + p4 )
pp1/pp2
dev.off()


# cell proportion distribution 
pdf("./04_Macrophages/Macrophages_subtype_proportion.pdf",width=3,height=5)
df <- as.data.frame(Macrophages@meta.data)
df_ct <- df %>%
    group_by(orig.ident, subtype) %>%
    summarise(counts = n()) %>%
    mutate(cell_proportion = counts / sum(counts))

p <- ggplot(df_ct, aes(orig.ident, cell_proportion, fill=subtype)) + 
    geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = my47colors[-1:-6]) +
    theme_cowplot() +
    xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
p
dev.off()


pdf("./04_Macrophages/Macrophages_subtype_M.pulmonis_violinplot.pdf",width=6,height=8)

p1 <- VlnPlot(subset(Macrophages,atac_S_M.pulmonis_ReadsPercentage>0), "atac_S_M.pulmonis_ReadsPercentage",add.noise=FALSE, pt.size = 0,col=my47colors[-1:-2],group.by = "subtype") +
geom_boxplot(width=0.1,outlier.size=0.01, position = position_dodge(0.9))

p2 <- VlnPlot(Macrophages, "atac_S_M.pulmonis_ReadsPercentage", pt.size = 0.01,col=my47colors[-1:-2],group.by = "subtype") +
geom_boxplot(width=0.1,outlier.size=0.01, position = position_dodge(0.9))
p1/p2
dev.off()




