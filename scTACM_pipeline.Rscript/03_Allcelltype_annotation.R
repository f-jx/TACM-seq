## create Seurat by gex matrix 
.libPaths("/data/R02/nieyg/ori/biosoft/conda/envs/r4-base/lib/R/library")

#library(BSgenome)
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

OSN_mm<-readRDS("./02_All_celltype/WNN_OSN_mm_integrated_all_celltype.rds")


########## 这里提取各主成分，并绘制不同Cluster的距离tree图，用来协助细胞分类，需要学习
###### marker gene annotation ######
# RNA tree 

Idents(OSN_mm)<- OSN_mm$seurat_clusters
object<- OSN_mm
embeddings <- Embeddings(object = object, reduction = "pca")[,1:50]
data.dims <- lapply(X = levels(x = object), FUN = function(x) {
    cells <- WhichCells(object = object, idents = x)
    if (length(x = cells) == 1) {
        cells <- c(cells, cells)
    }
    temp <- colMeans(x = embeddings[cells, ])
})
data.dims <- do.call(what = "cbind", args = data.dims)
colnames(x = data.dims) <- levels(x = object)
library(lsa)
cosine_dist <- as.dist(1-cosine(data.dims))
data.tree <- ape::as.phylo(x = hclust(d = cosine_dist))
library(ggtree);
pdf("./02_All_celltype/WNN_cluster_RNA_tree.pdf",width=5,height=5)
ggtree(data.tree,layout = "circular") + geom_tiplab()+ geom_treescale()
dev.off()


# ATAC tree 

embeddings <- Embeddings(object = object, reduction = "integrated_LSI")[,2:50]
data.dims <- lapply(X = levels(x = object), FUN = function(x) {
    cells <- WhichCells(object = object, idents = x)
    if (length(x = cells) == 1) {
        cells <- c(cells, cells)
    }
    temp <- colMeans(x = embeddings[cells, ])
})
data.dims <- do.call(what = "cbind", args = data.dims)
colnames(x = data.dims) <- levels(x = object)
library(lsa)
cosine_dist <- as.dist(1-cosine(data.dims))
data.tree <- ape::as.phylo(x = hclust(d = cosine_dist))
library(ggtree);
pdf("./02_All_celltype/WNN_cluster_ATAC_tree.pdf",width=5,height=5)
ggtree(data.tree,layout = "circular") + geom_tiplab()+ geom_treescale()
dev.off()


#####Annotate cells by RNA assay################
DefaultAssay(OSN_mm) <- "RNA"
Idents(OSN_mm)<-OSN_mm$seurat_clusters
# visualize the activities of canonical marker genes
markers <- c("Syt1", #神经元
              "Omp",# Mature ORNs
              "Nqo1","Ncam2","Cd36",
              "Gap43","Gng8",#Immature ORNs
              "Sox11","Neurod1","Neurog1",#INPs,immediate neuronal precursors 
              "Ascl1","Kit" ,#GBCs
              "Cebpd","Krt5","Trp63","Krt14",#HBCs
              "Sox2","Ermn","Cyp2g1","Cyp1a2",#支持细胞  Sustentacular cells
              "Atp1a2","Fabp7", #Ensheathing glia
              "Sox9","Sox10",#Bowman's gland
              "Pebp1","Calb2", #球周细胞Periglomerular cells
              "Ascl3","Cftr",#Microvillar cells
              "Krt18","Trpm5", #Brush cells
              "Col1a1","Bglap",#Osteogenic cells 成骨细胞
              "Eng","Sox17",#Pericytes
              "Cd37","Cd79a",#B cells
              "S100a9","S100a8",#中性粒细胞 Neutrophils
              "Hmgb2","Top2a",#Late activated neural stem cells
              "Lyz2","S100a4", #Monocytes
              "Hbb-bs","Hbb-bt",#红细胞,Erythrocytes
              "Mcpt8","Ccl4","Ptprc","Itga2","Fcer1a", #Basophils
              "C1qa","Ms4a7"#Macrophages
)

# RNA 
DefaultAssay(OSN_mm) <- 'RNA'
pdf("./02_All_celltype/OSN_mm_cluster-RNA_annotation-all_celltype.pdf",width=12,height=8)
p<-DotPlot(OSN_mm, features = markers,dot.scale = 3)+theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
p
dev.off()

# ATAC
library(BSgenome.Mmusculus.UCSC.mm10)
# ORN recall peak for subcluster 
DefaultAssay(OSN_mm)<-"ATAC"
peak<-CallPeaks(
       OSN_mm,
       group.by = "seurat_clusters",
       macs2.path = "/public/home/nieyg/biosoft/conda/bin/macs2",
       broad = FALSE,
       format = "BED",
       fragment.tempdir = tempdir(),
       effective.genome.size = 1.87e+09, #mm10 size
       outdir="./02_All_celltype/",
       combine.peaks=TRUE
)
macs2_counts <- FeatureMatrix(
     fragments = Fragments(OSN_mm),
     features = peak,
     cells = colnames(OSN_mm)
     )     
OSN_mm[["peaks_All_cluster"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  #fragments = fragpath[[i]],
  fragments = Fragments(OSN_mm),
  annotation = Annotation(OSN_mm)
)

########## 这里从染色质可及性数据计算推断基因转录活性，将scATAC矩阵转换为伪scRNA的数据（有染色质可及性推断的mRNA表达，非真实），计算各Cluster的Marker gene表达，辅助判断细胞类型，需要学习
DefaultAssay(OSN_mm) <- 'peaks_All_cluster'
gene.activities <- GeneActivity(OSN_mm)
OSN_mm[['ACTIVITY']] <- CreateAssayObject(counts = gene.activities)
OSN_mm <- NormalizeData(
  object = OSN_mm,
  assay = 'ACTIVITY',
  normalization.method = 'LogNormalize',
  scale.factor = median(OSN_mm$nCount_ACTIVITY)
)

DefaultAssay(OSN_mm) <- 'ACTIVITY'
pdf("./02_All_celltype/OSN_mm_cluster-ACTIVITY_annotation-all_celltype.pdf",width=12,height=8)
p<-DotPlot(OSN_mm, features = markers,dot.scale = 3)+theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
p
dev.off()


saveRDS(OSN_mm,"./02_All_celltype/WNN_OSN_integrated_all_celltype_recallpeaks.rds")

# Trackplot for annotation 
##Track for Marker genes promoters
Idents(OSN_mm)<-OSN_mm$seurat_clusters
library(BSgenome.Mmusculus.UCSC.mm10)
DefaultAssay(OSN_mm) <- "peaks_All_cluster"
# first compute the GC content for each peak
# conda r4-base环境有问题，用 r43计算
conda activate r43
.libPaths("/data/R02/nieyg/ori/biosoft/conda/envs/r43/lib/R/library")
library(Matrix)
library(Signac)
library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(sctransform)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(BSgenome.Mmusculus.UCSC.mm10)
set.seed(1234)
setwd("/data/R02/nieyg/project/OSN_mm/joint")
OSN_mm<-readRDS("./02_All_celltype/WNN_OSN_integrated_all_celltype_recallpeaks.rds")
DefaultAssay(OSN_mm) <- 'peaks_All_cluster'
########## 这一步是为了将染色质可及性与gene表达关联链接起来，从而绘制每个Cluster每个Marker基因的peak+gene表达图，辅助判断细胞类型，需要学习
OSN_mm <- RegionStats(OSN_mm, genome = BSgenome.Mmusculus.UCSC.mm10)
saveRDS(OSN_mm,"./02_All_celltype/WNN_OSN_integrated_all_celltype_recallpeaks_RegionStats.rds")

OSN_mm_test<- readRDS("./02_All_celltype/WNN_OSN_integrated_all_celltype_recallpeaks_RegionStats.rds")
#Annotation(OSN_mm)$tx_id <- Annotation(honeybee)$gene_name 

# link peaks to genes
OSN_mm <- LinkPeaks(
  object = OSN_mm,
  peak.assay = "peaks_All_cluster",
  expression.assay = "RNA",
  genes.use = markers
)
######Visulize track and RNA exp######
idents.plot <- Idents(OSN_mm)
markers <- c(
              "Cd37","Cd79a",#B cells
              "S100a9","S100a8",#中性粒细胞 Neutrophils
              "Hmgb2","Top2a",#Late activated neural stem cells
              "Lyz2","S100a4", #Monocytes
              "Hbb-bs","Hbb-bt",#红细胞,Erythrocytes
              "Mcpt8","Ccl4","Ptprc","Itga2","Fcer1a", #Basophils
              "C1qa","Ms4a7"#Macrophages
)

pdf("./02_All_celltype/Marker_gene-peaktrack-RNAexp-WNN.pdf",height=8,width=8)
for(i in markers){
  print(i)
  p1 <- CoveragePlot(
  object = OSN_mm,
  region = i,
  features = i,
  expression.assay = "RNA",
  idents = idents.plot,
  extend.upstream = 500,
  annotation=TRUE,
  extend.downstream = 500
)
print(p1)}
dev.off()


Idents(OSN_mm)<-OSN_mm$seurat_clusters
#####further annotation########
OSN_mm <- RenameIdents(
  object = OSN_mm,
  '0' = 'Neutrophils',
  '1' = 'Mature OSNs',
  '2' = 'B cells',
  '3' = 'Neutrophils',
  '4' = 'Basophils1',
  '5' = 'Macrophages',
  '6' = 'Neutrophils',
  '7' = 'Immature ORNs',
  '8' = 'HBCs',
  '9' = 'Sustentacular cells',
  '10' = "Bowman's gland",
  '11' = 'Sustentacular cells',
  '12' = 'Sustentacular cells',
  '13' = 'HBCs',
  '14' = 'Basophils2',
  '15' = 'Macrophages',
  '16' = 'B cells',
  '17' = 'Erythrocytes',
  '18' = 'Periglomerular cells',
  '19' = 'Monocytes',
  '20' = 'B cells',
  '21' = 'Microvillar cells',
  '22' = 'Sustentacular cells',
  '23' = 'Brush cells',
  '24' = "Bowman's gland",
  '25' = 'Mature OSNs',
  '26' = 'Mature OSNs'
  )
OSN_mm@meta.data$cell_type<-Idents(OSN_mm)
table(OSN_mm$cell_type,OSN_mm$orig.ident)

OSN_mm$cell_type<-factor(OSN_mm$cell_type,
  levels=c("HBCs","Immature ORNs","Mature OSNs",
    "Sustentacular cells","Microvillar cells","Bowman's gland","Periglomerular cells",
    "Brush cells","Basophils1","Basophils2","B cells","Erythrocytes","Monocytes","Macrophages","Neutrophils"
    ))
Idents(OSN_mm)<-OSN_mm$cell_type;

# Dimplot with annotation
pdf("./02_All_celltype/OSN_annotation_allcelltype_UMAP.pdf",width=6,height=5)
DimPlot(OSN_mm, label = T, repel = TRUE, cols=my47colors, reduction = "wnn.umap",group.by = "cell_type")
DimPlot(OSN_mm, label = F, repel = TRUE, cols=my47colors, reduction = "wnn.umap",group.by = "cell_type")+ ggtitle("")
dev.off()

# Save rds have annotation information 
DefaultAssay(OSN_mm) <- "RNA"
saveRDS(OSN_mm,"./02_All_celltype/WNN_OSN_integrated_all_celltype.rds")


# Dotplot with annotation
DefaultAssay(OSN_mm) <- 'RNA'
markers <- c("Syt1", #神经元
              "Omp",# Mature ORNs
              "Nqo1","Ncam2","Cd36",
              "Gap43","Gng8",#Immature ORNs
              "Cebpd","Krt5","Trp63","Krt14",#HBCs
              "Sox2","Ermn","Cyp2g1","Cyp1a2",#支持细胞  Sustentacular cells
              "Sox9","Sox10",#Bowman's gland
              "Pebp1","Calb2", #球周细胞Periglomerular cells
              "Ascl3","Cftr",#Microvillar cells
              "Krt18","Trpm5", #Brush cells
              "Cd37","Cd79a",#B cells
              "S100a9","S100a8",#中性粒细胞 Neutrophils
              "Lyz2","S100a4", #Monocytes
              "Hbb-bs","Hbb-bt",#红细胞,Erythrocytes
              "Mcpt8","Ccl4","Ptprc","Itga2","Fcer1a", #Basophils
              "C1qa","Ms4a7"#Macrophages
)

pdf("./02_All_celltype/OSN_mm_celltype-RNA_annotation.pdf",width=8,height=4)
p<-DotPlot(OSN_mm, features = markers,dot.scale = 3)+theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
p
dev.off()

pdf("./02_All_celltype/OSN_mm_celltype-ACTIVITY_annotation.pdf",width=8,height=4)
DefaultAssay(OSN_mm)<- "ACTIVITY"
p<-DotPlot(OSN_mm, features = markers,dot.scale = 3)+theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
p
dev.off()

# Trackplot with annotation

# conda r4-base环境有问题，用 r43计算
conda activate r43
.libPaths("/data/R02/nieyg/ori/biosoft/conda/envs/r43/lib/R/library")
library(Matrix)
library(Signac)
library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(sctransform)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(BSgenome.Mmusculus.UCSC.mm10)
set.seed(1234)
setwd("/data/R02/nieyg/project/OSN_mm/joint")
OSN_mm<-readRDS("./02_All_celltype/WNN_OSN_integrated_all_celltype.rds")
DefaultAssay(OSN_mm) <- 'peaks_All_cluster'

# link peaks to genes
OSN_mm <- LinkPeaks(
  object = OSN_mm,
  peak.assay = "peaks_All_cluster",
  expression.assay = "RNA",
  genes.use = markers
)
######Visulize track and RNA exp######
idents.plot <- Idents(OSN_mm)
markers <- c("Syt1", #神经元
              "Omp",# Mature ORNs
              "Nqo1","Ncam2","Cd36",
              "Gap43","Gng8",#Immature ORNs
              "Trp63","Krt14",#HBCs
              "Sox2","Ermn","Cyp2g1","Cyp1a2",#支持细胞  Sustentacular cells
              "Sox9","Sox10",#Bowman's gland
              "Pebp1","Calb2", #球周细胞Periglomerular cells
              "Ascl3","Cftr",#Microvillar cells
              "Krt18","Trpm5", #Brush cells
              "Cd37","Cd79a",#B cells
              "S100a9","S100a8",#中性粒细胞 Neutrophils
              "Lyz2","S100a4", #Monocytes
              "Hbb-bs","Hbb-bt",#红细胞,Erythrocytes
              "Mcpt8","Ccl4","Ptprc","Itga2","Fcer1a", #Basophils
              "C1qa","Ms4a7"#Macrophages
)

pdf("./02_All_celltype/All_Marker_gene-peaktrack_celltype.pdf",height=8,width=8)
for(i in markers){
  print(i)
  p1 <- CoveragePlot(
  object = OSN_mm,
  region = i,
  features = i,
  expression.assay = "RNA",
  idents = idents.plot,
  extend.upstream = 500,
  annotation=TRUE,
  extend.downstream = 500
)
print(p1)}
dev.off()


