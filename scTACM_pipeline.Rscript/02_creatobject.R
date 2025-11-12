
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
#library(DropletUtils)

set.seed(1234)

# create a Seurat object containing the filtered RNA and ATAC data;
# 1. RNA filtered matrix 
positive_counts <- Read10X_h5("/md01/nieyg/project/OSN_mm/data/cellranger/Positive/outs/filtered_feature_bc_matrix.h5")
positive_fragpath <- "/md01/nieyg/project/OSN_mm/data/cellranger/Positive/outs/atac_fragments_filtered.tsv.gz"
########## 这个是scRNA进一步QC的产物，需要QC源代码
positive_filtered_RNA_data<-read.table("/md01/nieyg/project/OSN_mm/data/cellranger/Positive/outs/counts_no_double_umi_001.tsv.gz")
positive_filtered_RNA_counts = spread(positive_filtered_RNA_data, V3, V2)
positive_filtered_RNA_counts[is.na(positive_filtered_RNA_counts)] <- 0
rownames(positive_filtered_RNA_counts)<-positive_filtered_RNA_counts$V1
positive_filtered_RNA_counts<-positive_filtered_RNA_counts[,-1]
positive_RNA_barcode<-paste(colnames(positive_filtered_RNA_counts),"-1",sep="")
colnames(positive_filtered_RNA_counts)<-positive_RNA_barcode

negative_counts <- Read10X_h5("/md01/nieyg/project/OSN_mm/data/cellranger/Negative/outs/filtered_feature_bc_matrix.h5")
negative_fragpath <- "/md01/nieyg/project/OSN_mm/data/cellranger/Negative/outs/atac_fragments_filtered.tsv.gz"
negative_filtered_RNA_data<-read.table("/md01/nieyg/project/OSN_mm/data/cellranger/Negative/outs/counts_no_double_umi_001.tsv.gz")
negative_filtered_RNA_counts = spread(negative_filtered_RNA_data, V3, V2)
negative_filtered_RNA_counts[is.na(negative_filtered_RNA_counts)] <- 0
rownames(negative_filtered_RNA_counts)<-negative_filtered_RNA_counts$V1
negative_filtered_RNA_counts<-negative_filtered_RNA_counts[,-1]
negative_RNA_barcode<-paste(colnames(negative_filtered_RNA_counts),"-1",sep="")
colnames(negative_filtered_RNA_counts)<-negative_RNA_barcode

# 2. ATAC filtered matrix 

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))
########## 读取peaks.bed文件并merge，作为ATAC矩阵的行
region_positive <- read.table(
  file = "/md01/nieyg/project/OSN_mm/data/cellranger/Positive/outs/atac_peaks.bed",
  col.names = c("chr", "start", "end")
)
region_negative <- read.table(
  file = "/md01/nieyg/project/OSN_mm/data/cellranger/Negative/outs/atac_peaks.bed",
  col.names = c("chr", "start", "end")
)
positive_gr <- makeGRangesFromDataFrame(region_positive)
negative_gr <- makeGRangesFromDataFrame(region_negative)
combined.region <- GenomicRanges::reduce(x = c(positive_gr, negative_gr))
# Filter out bad peaks based on length
########## 这里对merge后的peak的宽度做了筛选
widths <- width(combined.region)
combined.region <- combined.region[widths  < 10000 & widths > 20]
combined.region

positive_filtered_frag <- CreateFragmentObject(
  path = positive_fragpath,
  cells = positive_RNA_barcode,
  validate.fragments = FALSE
)
positive_filtered_chrom_matrix<-FeatureMatrix(
  fragments = positive_filtered_frag,
  features = combined.region
)
positive_filtered_chrom_assay <- CreateChromatinAssay(
  counts = positive_filtered_chrom_matrix,
  sep = c(":", "-"),
  fragments = positive_filtered_frag,
  annotation=annotation
)

negative_filtered_frag <- CreateFragmentObject(
  path = negative_fragpath,
  cells = negative_RNA_barcode,
  validate.fragments = FALSE
)
negative_filtered_chrom_matrix<-FeatureMatrix(
  fragments = negative_filtered_frag,
  features = combined.region
)
negative_filtered_chrom_assay <- CreateChromatinAssay(
  counts = negative_filtered_chrom_matrix,
  sep = c(":", "-"),
  fragments = negative_filtered_frag,
  annotation=annotation
)

# 3. create object 

positive <-  CreateSeuratObject(counts = positive_filtered_RNA_counts, min.cells=3, assay = "RNA")
positive[["ATAC"]]<- positive_filtered_chrom_assay

negative <-  CreateSeuratObject(counts = negative_filtered_RNA_counts, min.cells=3, assay = "RNA")
negative[["ATAC"]]<- negative_filtered_chrom_assay

objList<- list(positive,negative)
# calculate the score of NS and TSS
for (i in seq_len(length(objList))) {
  DefaultAssay(objList[[i]]) <- "ATAC"
    objList[[i]] <- NucleosomeSignal(objList[[i]])
    objList[[i]] <- TSSEnrichment(objList[[i]],fast=FALSE)
    }
library(ggplot2)
pdf("./01_QC/DensityScatter_QC.pdf")
########## 这里的nCount_ATAC是在创建Seurat obj的过程中自动创建的。
DensityScatter(objList[[1]], x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
DensityScatter(objList[[2]], x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
dev.off()

# plot TSS and fragrament distribution plot 
  pdf("./01_QC/TSS_distribution_Positive.pdf")
  objList[[1]]$high.tss<-ifelse(objList[[1]]$TSS.enrichment > 4, 'High', 'Low')
  TSS<-TSSPlot(objList[[1]], group.by = 'high.tss') + NoLegend()+ labs(title = "NE")
  objList[[1]]$nucleosome_group <- ifelse(objList[[1]]$nucleosome_signal > 2, 'NS > 2', 'NS < 2')
   print(TSS);
  dev.off();
  pdf("./01_QC/TSS_distribution_Negative.pdf")
  objList[[2]]$high.tss<-ifelse(objList[[2]]$TSS.enrichment > 4, 'High', 'Low')
  TSS<-TSSPlot(objList[[2]], group.by = 'high.tss') + NoLegend()+ labs(title = "Nurse")
  objList[[2]]$nucleosome_group <- ifelse(objList[[2]]$nucleosome_signal > 2, 'NS > 2', 'NS < 2')
  print(TSS);
  dev.off();
pdf("./01_QC/FragmentHistogram_QC.pdf")
FragmentHistogram(objList[[1]], group.by = 'nucleosome_group', region = "chr1-1-20000000")
FragmentHistogram(objList[[2]], group.by = 'nucleosome_group', region = "chr1-1-20000000")
dev.off()

sample<- c("Positive","Negative")
for (i in seq_len(length(objList))) {
  # plot QC plot 
  pdf(file = paste("./01_QC/3b_QC_before_",sample[i],".pdf", sep = ""),width=10,height=6)
  qc<-VlnPlot(object = objList[[i]],
            features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal","percent.mt","nFeature_RNA"),
            ncol = 5,
            pt.size = 0.01
          )
  densityRNA<-plot(density(objList[[i]]@meta.data$nCount_RNA),xlim=c(0,1000))
  densityRNA_feature<-plot(density(objList[[i]]@meta.data$nFeature_RNA),xlim=c(0,500))
  densityATAC<-plot(density(objList[[i]]@meta.data$nCount_ATAC),xlim=c(0,10000))
  print(qc)
  print(densityRNA)
  print(densityRNA_feature)
  print(densityATAC)
  dev.off();
   }


   # filter out low quality cells
# To remove doublets,select different cutoff#####
  objList2<-list()
for (i in 1:2){
  objList2[[i]] <- subset(objList[[i]],
     subset= nCount_ATAC < 10000 & nCount_ATAC > 200 &
      nCount_RNA < 5000 & nCount_RNA > 50 &
      nFeature_RNA < 1000 & nFeature_RNA > 20 &
      TSS.enrichment >= 4 &
      nucleosome_signal < 2 
      )
}


  objList3<-list()
for (i in 1:2){
  objList3[[i]] <- subset(objList[[i]],
     subset= nCount_ATAC < 10000 & nCount_ATAC > 200 &
      TSS.enrichment >= 4 &
      nucleosome_signal < 2 
      )
}



objList<-objList2

# Peak Calling 
# quantify counts in each peak
fragpath<-list(positive_fragpath,negative_fragpath)
for(i in seq_len(length(objList))){
  ########## 重新Callpeak并计算fragments矩阵
	# call peaks using MACS2
    peaks <- CallPeaks(objList[[i]], macs2.path = "/public/home/nieyg/biosoft/conda/bin/macs2")
    macs2_counts <- FeatureMatrix(
      fragments = Fragments(objList[[i]]),
      features = peaks,
      cells = colnames(objList[[i]])
    )     
    # create a new assay using the MACS2 peak set and add it to the Seurat object
    objList[[i]][["peaks"]] <- CreateChromatinAssay(
      counts = macs2_counts,
      fragments = fragpath[[i]],
      annotation = annotation
    )
}

#######integrate RNA and ATAC#####################
# Simply merge Seurat objects
########## 简单地merge只会简单地merge peak然后按原来各自的peak表达量来填充新矩阵，并不会按新peak重新计算表达量。
merged_obj <- merge(x=objList[[1]],y=objList[[2]],add.cell.ids = c("Positive","Negative"),project = "OSN_mm")
Idents(merged_obj) <- gsub("_.*", "", colnames(merged_obj))
merged_obj$orig.ident<-Idents(merged_obj)
OSN_mm.list <- SplitObject(merged_obj, split.by = "orig.ident")


# for RNA 
for (i in 1:length(OSN_mm.list)) {
	DefaultAssay(OSN_mm.list[[i]]) <- "RNA"
  OSN_mm.list[[i]] <- SCTransform(OSN_mm.list[[i]], verbose = FALSE)
}
for (i in seq_len(length(OSN_mm.list))) {
  DefaultAssay(OSN_mm.list[[i]]) <- "SCT"
}
########## RNA的integrate去批次方案，这里使用SCT+传统的寻找共同锚点的方式，没有用IntegrateLayers的方式。
OSN_mm.features <- SelectIntegrationFeatures(object.list = OSN_mm.list, nfeatures = 1000)
OSN_mm.list <- PrepSCTIntegration(object.list = OSN_mm.list, anchor.features = OSN_mm.features)
#integrate RNA using rpca
OSN_mm_list <- lapply(
  X = OSN_mm.list,
  FUN = RunPCA,
  features = OSN_mm.features,
  verbose = FALSE
)

integration_anchors_RNA <- FindIntegrationAnchors(
  object.list = OSN_mm_list,
  normalization.method = "SCT",
  anchor.features = OSN_mm.features,
  dims = 1:30,
  reduction = "rpca",
  k.anchor = 20,
)


OSN_mm_RNA <- IntegrateData(
  anchorset = integration_anchors_RNA,
  normalization.method = "SCT",
  new.assay.name = "integrated_RNA",
  dims = 1:30,
  features.to.integrate =OSN_mm.features
)
 
 # ATAC 
########## ATAC的integrate去批次方案，只能用传统的寻找共同锚点的方式
#run LSI on new seurat object with integrated RNA assay
# preprocess each object
for (i in 1:length(OSN_mm.list)){
  DefaultAssay(OSN_mm.list[[i]]) <- "ATAC"
  OSN_mm.list[[i]] <- RunTFIDF(OSN_mm.list[[i]])
  OSN_mm.list[[i]] <- FindTopFeatures(OSN_mm.list[[i]], min.cutoff = 20)
  OSN_mm.list[[i]] <- RunSVD(object=OSN_mm.list[[i]])
}
integration_anchors_ATAC <- FindIntegrationAnchors(
  object.list = list(OSN_mm.list[[1]], OSN_mm.list[[2]]),
  anchor.features = rownames(OSN_mm.list[[2]]),
  reduction = "rlsi",
  dims = 2:30
)

########## 这里ATAC简单merge的结果没有用common.peak.region重新计算表达矩阵，可能不合适
# process merged object
OSN_mm<- merged_obj
DefaultAssay(OSN_mm) <- "ATAC"
OSN_mm <- RunTFIDF(OSN_mm)
OSN_mm <- FindTopFeatures(OSN_mm, min.cutoff = 50)
OSN_mm <- RunSVD(object=OSN_mm)

OSN_mm_atac <- IntegrateEmbeddings(
  anchorset = integration_anchors_ATAC,
  new.reduction.name = "integrated_LSI",
  reductions = OSN_mm@reductions$lsi,
  dims.to.integrate = 1:30
)

########## 这里的做法我们一致，直接把原来的reduction和assays复制过来即可
#copy integrated LSI from duplicate seurat object to original object
OSN_mm@reductions$integrated_LSI <- OSN_mm_atac@reductions$integrated_LSI
OSN_mm@assays$integrated_RNA <- OSN_mm_RNA@assays$integrated_RNA

#####done integrate ATAC and RNA ################

# RNA analysis
DefaultAssay(OSN_mm) <- "integrated_RNA"
########## 这里RunPCA之前可能应该先运行ScaleData
OSN_mm <- RunPCA(OSN_mm) 
pdf("./02_All_celltype/ElbowPlot_RNA.pdf")
ElbowPlot(OSN_mm,ndims=50,reduction='pca')
dev.off();

OSN_mm <- RunUMAP(OSN_mm, dims = 1:30, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(OSN_mm) <- "ATAC"
pdf("./02_All_celltype/LSI-depth-correalation.pdf")
########## 这里绘图展示PCA各主成分和测序深度的相关性，剔除相关性高的主成分。
DepthCor(OSN_mm)
dev.off();
OSN_mm <- RunUMAP(OSN_mm, reduction = 'integrated_LSI', dims = 2:40, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

# build a joint neighbor graph using both assays
########## 加权合并ATAC和RNA，绘制WNN的UMAP图并分群
OSN_mm <- FindMultiModalNeighbors(
  object = OSN_mm,
  reduction.list = list("pca", "integrated_LSI"), 
  dims.list = list(1:30, 2:40),
  #modality.weight.name = "RNA.weight",
  verbose = TRUE
)
OSN_mm <- RunUMAP(OSN_mm, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
OSN_mm <- FindClusters(OSN_mm, graph.name = "wsnn", resolution =0.5, algorithm = 3, verbose = FALSE)

###reorder the level of sample#####
Idents(OSN_mm)<-OSN_mm$orig.ident
OSN_mm$orig.ident<-factor(OSN_mm$orig.ident,levels=c("Negative","Positive"))
########## 这里的颜色方案可以借用一下，很好看
my47colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
         '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
         '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
         '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', '#116530', '#678983',
         '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76')
pdf("./02_All_celltype/OSN_mm_cluster_WNN.pdf",width=15,height=5)
###cluster
p1 <- DimPlot(OSN_mm, cols=my47colors, reduction = "umap.rna", group.by = "seurat_clusters", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(OSN_mm, cols=my47colors, reduction = "umap.atac",group.by = "seurat_clusters", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(OSN_mm, cols=my47colors, reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("WNNUMAP")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
###sample
p1 <- DimPlot(OSN_mm, cols=my47colors, reduction = "umap.rna", group.by = "orig.ident", shuffle=TRUE,label = F, label.size = 5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(OSN_mm, cols=my47colors, reduction = "umap.atac",group.by = "orig.ident", shuffle=TRUE,label = F, label.size = 5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(OSN_mm, cols=my47colors, reduction = "wnn.umap", group.by = "orig.ident", shuffle=TRUE,label = F, label.size = 5, repel = TRUE) + ggtitle("WNNUMAP")
p1 + p2 + p3 & theme(plot.title = element_text(hjust = 0.5))#& NoLegend()
dev.off()


#add M.pulmonisReads in object 

M.pulmonisReads<- read.csv("/md01/nieyg/project/OSN_mm/joint/M.pulmonisReads.csv")
M.pulmonisReads$barcode <- gsub(
  pattern = "jointF1met11-19d", 
  replacement = "Positive",        
  x = M.pulmonisReads$barcode      
)
M.pulmonisReads$barcode <- gsub(
  pattern = "jointF1met14d100G2", 
  replacement = "Negative",        
  x = M.pulmonisReads$barcode      
)

head(M.pulmonisReads)
Negative_data2<- data.frame(barcode=colnames(OSN_mm.list[[2]]),atac_S_M.pulmonisReads=0,atac_S_M.pulmonis_ReadsPercentage=0,rna_S_M.pulmonisReads=0,rna_S_M.pulmonis_ReadsPercentage=0)
merge<- rbind(M.pulmonisReads,Negative_data2)
setdiff(colnames(OSN_mm),merge$barcode)
OSN_mm$atac_S_M.pulmonisReads<- merge[match(rownames(OSN_mm@meta.data),merge$barcode),]$atac_S_M.pulmonisReads
OSN_mm$atac_S_M.pulmonis_ReadsPercentage<- merge[match(rownames(OSN_mm@meta.data),merge$barcode),]$atac_S_M.pulmonis_ReadsPercentage
OSN_mm$rna_S_M.pulmonisReads<- merge[match(rownames(OSN_mm@meta.data),merge$barcode),]$rna_S_M.pulmonisReads
OSN_mm$rna_S_M.pulmonis_ReadsPercentage<- merge[match(rownames(OSN_mm@meta.data),merge$barcode),]$rna_S_M.pulmonis_ReadsPercentage


pdf('./02_All_celltype/atac_S_M.pulmonis_ReadsPercentage_FeaturePlot_WNN.pdf', width=15, height=5)
p1<- FeaturePlot(OSN_mm, reduction = 'umap.rna',features = "atac_S_M.pulmonis_ReadsPercentage", ncol = 1)+ ggtitle("RNA")
p2<- FeaturePlot(OSN_mm, reduction = 'umap.atac',features = "atac_S_M.pulmonis_ReadsPercentage", ncol = 1)+ ggtitle("ATAC")
p3<- FeaturePlot(OSN_mm, reduction = 'wnn.umap', features = "atac_S_M.pulmonis_ReadsPercentage", ncol = 1)+ ggtitle("WNN")
p1 | p2 | p3 
dev.off()

saveRDS(OSN_mm,"./02_All_celltype/WNN_OSN_mm_integrated_all_celltype.rds")












