
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

OSN_mm<- readRDS("./02_All_celltype/WNN_honeybee_integrated_all_celltype.rds")
########## 这里的配色方案好用
scRNA_newpalette <- c(brewer.pal(9,"YlGn")[c(9,7,5,4,2)],brewer.pal(12,"Set3")[8],brewer.pal(9,"RdPu")[c(2,4,5,7)],brewer.pal(9,"Blues")[c(2,6,8)],brewer.pal(9,"BuPu")[5],brewer.pal(9,"BuPu")[7])
scATAC_newpalette <- darken(scRNA_newpalette, 0.2)

# UMAP with annotation or Sample
cell_number_df <- integrated_scATAC@meta.data[,"cell_type",drop=FALSE] %>%
  group_by(cell_type) %>%
  summarise(counts=n())
cell_number_df$log10_counts <- log10(cell_number_df$counts)
cell_number_df$cell_type <- factor(cell_number_df$cell_type,levels=rev(c("HBCs","GBCs","INPs","Immature OSNs","Mature OSNs","Sustentacular cells","Ensheathing glia","Bowman's gland","Periglomerular cells","Brush cells","B cells","Erythrocytes","Basophils","Neutrophils","Monocytes")))
cell_type_df <- cell_number_df[,"cell_type",drop=FALSE]
cell_type_df$x <- 1
pdf("/md01/shipy3/Projects/mouse_ORs/output/scATAC/two_sample_OE_add/integrated_OE_scATAC_UMAP_with_cell_types.pdf",width=9.5)
p1 <- DimPlot(integrated_scATAC, reduction = "umap",label=TRUE,cols=scATAC_newpalette,repel=TRUE,label.size=5)+
  labs(title="scATAC-seq")+
  #geom_segment(aes(x = -14.641, y = -16.31245, xend = -12.2, yend = -16.31245),arrow = arrow(length = unit(0.3, "cm")))+
  #geom_segment(aes(x = -14.641, y = -16.31245, xend = -14.641, yend = -13.31245),arrow = arrow(length = unit(0.3, "cm")))+
  #annotate(geom = "text", x = -13.5, y = -17, label = "UMAP_1", color = "black") +
  #annotate(geom = "text", x = -15.2, y = -15, label = "UMAP_2", color = "black",angle = 90) +
  theme(plot.title = element_text(hjust = 0.5,size=rel(1.5),face="bold"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background=element_rect(fill='transparent', color='black',linetype="solid",size=1.5),axis.ticks=element_blank(),axis.text=element_blank(),axis.title=element_text(size=rel(1)),legend.position="none",plot.margin=unit(c(0.2, 0, 0.2, 0.3), "cm"))
p2 <- ggplot(data=cell_number_df,aes(x=log10_counts,y=cell_type,fill=cell_type))+
  geom_bar(stat = "identity",width=0.8)+
  labs(x="Log10 (# cells)")+
  #scale_x_continuous(limits=c(0,11100),breaks = c(1000,11000))+
  #scale_y_discrete(position = "right")+
  scale_fill_manual(values=rev(scATAC_newpalette))+
  theme(plot.title = element_text(hjust = 0.5,size=rel(1.5),face="bold"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background=element_rect(fill='transparent', color='black',linetype="solid",size=1.5),axis.ticks.y=element_blank(),axis.text.y=element_blank(),axis.title.y=element_blank(),legend.position="none",plot.margin=unit(c(0.2, 0, 0.2, 0), "cm"),axis.text.x=element_text(size=rel(1),color="black"),axis.title.x=element_text(size=rel(1)))
p3 <- ggplot(data=cell_type_df,aes(x=x,y=cell_type))+
  geom_tile(aes(fill=cell_type),color="black")+
  scale_y_discrete(position = "right",expand=c(0,0))+
  scale_fill_manual(values=rev(scATAC_newpalette))+
  theme(plot.title = element_text(hjust = 0.5,size=rel(1.5),face="bold"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background=element_rect(fill='transparent', color='black',linetype="solid",size=1.5),axis.ticks.x=element_blank(),axis.text.x=element_blank(),axis.title.x=element_blank(),legend.position="none",plot.margin=unit(c(0.2, 0, 0.2, 0.2), "cm"),axis.text.y=element_text(size=rel(1.5),color="black"),axis.title.y=element_blank())
p1 + p2 + p3 + plot_layout(widths=c(8,1.5,0.2))
dev.off()

## marker gene track for show 
## Figure S4G - marker peaks of all cell types
markers <- c("Syt1", #神经元
              "Omp",# Mature ORNs
              "Nqo1","Ncam2",
              "Gap43","Gng8",#Immature ORNs
              "Sox11","Neurod1","Neurog1",#INPs,immediate neuronal precursors 
              "Ascl1","Kit" ,#GBCs
              "Krt5","Trp63","Krt14",#HBCs
              "Sox2","Ermn","Cyp2g1","Cyp1a2",#支持细胞  Sustentacular cells
              "Atp1a2","Fabp7", #Ensheathing glia
              "Sox9","Sox10",#Bowman's gland
              "Pebp1","Calb2", #球周细胞Periglomerular cells
              "Ascl3","Cftr",#Microvillar cells
              "Krt18","Trpm5", #Brush cells
              "Col1a1","Bglap",#Osteogenic cells
              "Eng","Sox17",#Pericytes
              "Cd37","Cd79a",#B cells
              "S100a9","S100a8",#中性粒细胞 Neutrophils
              "Hmgb2","Top2a",#Late activated neural stem cells
              "Lyz2","S100a4", #Monocytes
              "Hbb-bs","Hbb-bt",#红细胞,Erythrocytes
              "Mcpt8","Ccl4","Itga2", #Basophils
              "C1qa","Ms4a7"#Macrophages
)
pdf(str_c(out_dir,"integrated_OE_scATAC_cell_type_markers_CoveragePlot.pdf"))
for (i in 1:length(markers)){
  p <- CoveragePlot(
    object = integrated_scATAC,
    region = markers[i],
    extend.upstream = 1000,
    extend.downstream = 1000,
    peaks=FALSE
  )
  print(p)
}
dev.off()



regions <- c(
  "chr5-75574200-75575500",# Kit
  "chr2-79456000-79457300",# Neurod1
  "chr16-42277000-42279100",# Gap43
  "chr7-98144500-98146500",# Omp
  "chr2-58052364-58053864",# Ermn
  "chr1-172297064-172298564",# Atp1a2
  "chr15-79164000-79166000",# Sox10
  "chr8-110167000-110168500",# Calb2
  "chr7-143093500-143094250",# Trpm5
  "chr7-24897000-24897700",# Cd79a
  "chr7-103827500-103828300",# Hbb−bs
  #"chr14-56085000-56085500",# Mcpt8 
  "chr13-114844000-114848000", # Itga2
  "chr3-90695400-90696200",# S100a9
  "chr3-90603300-90604000"# S100a4
  )
names <- c("Kit","Neurod1","Gap43","Omp","Ermn","Atp1a2","Sox10","Calb2","Trpm5","Cd79a","Hbb-bs","Itga2","S100a9","S100a4")
cov_plot.ls <- lapply(1:length(regions),function(i){
  cov_plot <- CoveragePlot(object = integrated_scATAC,region=regions[i],peaks=FALSE,annotation = FALSE)
  cov_df <- cov_plot$data 
  cov_plot <- ggplot(data = cov_df,mapping = aes(x = position, y = coverage, fill = group)) +
    geom_area(stat = "identity") +
    geom_hline(yintercept = 0, size = 0.1) +
    facet_wrap(facets = ~group, strip.position = "left", ncol = 1) +
    ylim(c(0, signif(max(cov_df$coverage, na.rm = TRUE), digits = 2))) +
    theme_bw() + 
    labs(title=names[i]) +
    theme(plot.title=element_text(face = "bold",size = rel(1.5),hjust = 0.5),panel.background=element_rect(fill='transparent', color='black',linetype="solid"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "none",strip.text.y.left = element_blank(), strip.background = element_blank(),axis.text.y = element_blank(),panel.spacing.y = unit(x = 0, units = "line"),plot.margin=unit(c(0, 0,0,0), "cm"),axis.title.y = element_blank(),axis.ticks.y = element_blank(),axis.text.x = element_blank(),axis.title.x = element_blank(),axis.ticks=element_blank()) +
    scale_fill_manual(values = scATAC_newpalette)
  cov_plot
  })

gene_plot.ls <- lapply(1:length(regions),function(i){
  gene_plot <- AnnotationPlot(object = integrated_scATAC,region = regions[i])
  gene_plot <- gene_plot + theme(panel.background=element_rect(fill='transparent', color='black',linetype="solid"),axis.text.x = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank(),axis.title.y=element_blank(),plot.margin=unit(c(0.2,0.03,0,0), "cm"))
  gene_plot
  })

pdf(str_c(out_dir,"integrated_OE_scATAC_all_cell_type_markers_combined_CoveragePlot.pdf"),width=18,height=6.5)
((Krt5_cov_plot + Krt5_gene_plot + plot_layout(heights=c(7,0.6))) | 
(cov_plot.ls[[1]] + gene_plot.ls[[1]] + plot_layout(heights=c(7,0.6))) |
(cov_plot.ls[[2]] + gene_plot.ls[[2]] + plot_layout(heights=c(7,0.6))) |
(cov_plot.ls[[3]] + gene_plot.ls[[3]] + plot_layout(heights=c(7,0.6))) |
(cov_plot.ls[[4]] + gene_plot.ls[[4]] + plot_layout(heights=c(7,0.6))) |
(cov_plot.ls[[5]] + gene_plot.ls[[5]] + plot_layout(heights=c(7,0.6))) |
(cov_plot.ls[[6]] + gene_plot.ls[[6]] + plot_layout(heights=c(7,0.6))) |
(cov_plot.ls[[7]] + gene_plot.ls[[7]] + plot_layout(heights=c(7,0.6))) |
(cov_plot.ls[[8]] + gene_plot.ls[[8]] + plot_layout(heights=c(7,0.6))) |
(cov_plot.ls[[9]] + gene_plot.ls[[9]] + plot_layout(heights=c(7,0.6))) |
(cov_plot.ls[[10]] + gene_plot.ls[[10]] + plot_layout(heights=c(7,0.6))) |
(cov_plot.ls[[11]] + gene_plot.ls[[11]] + plot_layout(heights=c(7,0.6))) |
(cov_plot.ls[[12]] + gene_plot.ls[[12]] + plot_layout(heights=c(7,0.6))) |
(cov_plot.ls[[13]] + gene_plot.ls[[13]] + plot_layout(heights=c(7,0.6))) | 
(cov_plot.ls[[14]] + gene_plot.ls[[14]] + plot_layout(heights=c(7,0.6)))) +  plot_layout(widths=c(1.1,1,1,1,1,1,1,1,1,1,1,1,1,1,1))
dev.off()

## Figure S4B - UMAP split by sample
newpalette <- c("#F39B7FFF","#E64B35FF")
pdf(str_c(out_dir,"integrated_OE_scATAC_split_by_sample_UMAP.pdf"),width=7,height=12)
DimPlot(integrated_scATAC, reduction = "umap",label=FALSE,cols=newpalette,group.by="dataset",split.by="dataset",ncol=1)+labs(title="scATAC-seq")+theme(plot.title = element_text(hjust = 0.5,size=rel(1.5),face="bold"))
dev.off()





