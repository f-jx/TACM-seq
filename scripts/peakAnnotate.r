# ChIPseeker使用：
# https://www.jianshu.com/p/c76e83e6fa57
# https://www.jianshu.com/p/5fb041f09953
# 还有GREAT网站也可以对peak进行在线注释
# pdf(str_replace(f,'.xls','.coverage.split.pdf'),height=14,width=14)
#     print(ggarrange(
#         covplot(peak.up,title="Peaks UP"),
#         covplot(peak.dw,title="Peaks DOWN"),
#         ncol=1,nrow=2))
# dev.off()
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install(ask=FALSE)  # 更新 BiocManager 本身
# BiocManager::install(c("clusterProfiler","org.Hs.eg.db","org.Mm.eg.db","GO.db","KEGG.db","KEGGREST"),ask=FALSE) # 替换为你需要更新的包
args <- c('','~/project/refSeq/UCSC/mm39.ensGene.gtf','org.Mm.eg.db','mmu','enrich')
args <- commandArgs(trailingOnly = TRUE)
diff_peak.files <- args[[1]];if(diff_peak.files=='')quit(save='no')
gtf.file <- args[[2]]
orgdb <- args[[3]]
organism <- args[[4]]
enrich <- ifelse(length(args)>4,args[[5]],'')

library(tidyverse)
library(ggpubr)

get_peaks_flank_geneIds_transIds <- function(gr_peak,gr_promoters){
    overlaps <- findOverlaps(gr_peak,gr_promoters)
    hits <- data.frame(
        peak_index = queryHits(overlaps),
        tx_name = mcols(gr_promoters)$tx_name[subjectHits(overlaps)],
        gene_id = mcols(gr_promoters)$gene_id[subjectHits(overlaps)]
    ) %>%group_by(peak_index) %>%
        summarise(
            tx_name = paste(unique(tx_name), collapse = ";"),
            gene_id = paste(unique(gene_id), collapse = ";"))  # 居然可以这样！！！第一次见到
    return(hits)
}
# strsplit(geneID,';')[[1]]返回向量，去除点后缀，match返回每个元素对应在df$ENTREZID的位置（即行数）
# 再由行数定位到对应$SYMBOL和ENTREZID，并用';'连接起来，完成转换
get_flank_col<-function(x,df,keycol,targetcol)paste(df[[targetcol]][match(str_replace(strsplit(x,";")[[1]],'\\..*',''),df[[keycol]])],collapse=';')
del_NA_dup <- function(x) {
  elements <- strsplit(x, ";")[[1]]          # 分割字符串
  elements_clean <- elements[!elements %in% c(NA,"NA","")] # 删除NA
  if(length(elements_clean)==0)return('')
  elements_unique <- unique(elements_clean)   # 去重
  paste(elements_unique, collapse = ";")      # 合并
}
peakAnot <- function(name,peak,txdb,orgdb,flkdist=50000,promoter_upstream=2000,promoter_downstream=0){
    library(ChIPseeker)
    library(orgdb,character.only=TRUE)
    # peak注释：    https://mp.weixin.qq.com/s/vWTf59KDs1lp_sQXjEhI_g
    # 1功能区域annotation（promoter、UTR），由assignGenomicAnnotation控制，默认打开，做基因可变剪接的，关注peak所在基因组上的具体位置，选这个；
    # 2（nearest geneId、transcriptId），tssRegion设置的是要把转录起始位点上下游多长的区域标记为Promoter（annotation列），不会改变表格最终peak数量，也不会改变distance列的值；
    # 3（up、downstream flank注释），以peak为中心上下游搜索多少bp的距离将这区间的所有基因都考虑进来。
    # TxDb将bed区间注释成geneID，annoDb进一步通过geneID在末尾附加ENSEMBL、SYMBOL、GENENAME等几列信息
    if(file.exists(paste0(name,'.annotation.xls.gz'))){  # 如果文件存在，可以直接从文件中读取，加快速度
        print('########## 原注释文件存在，将从文件中读取信息以加快速度')
        anot <- annotatePeak(peak,assignGenomicAnnotation=TRUE,tssRegion=c(-promoter_upstream, promoter_downstream),TxDb=txdb,annoDb=orgdb)
        tb.anot <- anot@anno%>%as_tibble(); tb.file <- read_tsv(paste0(name,'.annotation.xls.gz'))
        tb <- left_join(tb.anot[,1:3],tb.file); selectcols <- !names(tb)%in%names(tb.anot)
        mcols(anot@anno) <- cbind(mcols(anot@anno),tb[,selectcols])
        return(anot)
    }
    anot <- annotatePeak(peak,assignGenomicAnnotation=TRUE,
        tssRegion = c(-promoter_upstream, promoter_downstream),TxDb =txdb,
        annoDb=orgdb,addFlankGeneInfo=TRUE,flankDistance=flkdist)
    # 统计transcript和gene的启动子区域上的peak
    gr_trans <- transcripts(txdb,columns=c('gene_id','tx_name'))
    # transcripts()获得的gene_id是GRangesList格式，需要转换成向量，使其和gr_gene一致
    gr_trans$gene_id <- gr_trans$gene_id@unlistData
    # mcols()函数可以获取GRanges对象的metadata，这里给GRange对象的metadata添加tx_name列，使其和gr_trans一致
    gr_gene <- genes(txdb); mcols(gr_gene)$tx_name <- 'NA'
    gr_promoters <- c(gr_trans,gr_gene)     # 将gene和trans两个GRange对象合并起来，准备计算启动子区域
    # 计算基因的启动子区域，genes()输出的结果，不论正负链start数值总小于end，正链基因TSS在start前，负链基因TSS在end后
    ranges(gr_promoters) <- IRanges(
        start=ifelse(strand(gr_promoters)=='+',start(gr_promoters)-promoter_upstream,end(gr_promoters)-promoter_downstream),
        end=ifelse(strand(gr_promoters)=='+',start(gr_promoters)+promoter_downstream,end(gr_promoters)+promoter_upstream)
    )
    hits <- get_peaks_flank_geneIds_transIds(anot@anno,gr_promoters)
    # mcols(anot@anno)等同于anot@anno@elementMetadata@listData，这里必须要新赋值，因为后面并不是都有内容，很多是注释不到的
    mcols(anot@anno)$promoter_flank_txIds <- rep(NA,anot@anno@elementMetadata@nrows)
    mcols(anot@anno)$promoter_flank_txIds[hits$peak_index] <- hits$tx_name;
    mcols(anot@anno)$promoter_flank_geneIds <- rep(NA,anot@anno@elementMetadata@nrows)
    mcols(anot@anno)$promoter_flank_geneIds[hits$peak_index] <- hits$gene_id;

    print('########## 将flank_geneIds转换成flank_symbols和flank_entrezIds')
    library(clusterProfiler)
    k <- keys(get(orgdb),keytype='ENSEMBL')
    df_ENSEMBL <- bitr(k,fromType='ENSEMBL',toType=c('SYMBOL','ENTREZID'),OrgDb=orgdb,drop = TRUE)
    k <- keys(get(orgdb),keytype='ENSEMBLTRANS')
    df_ENSEMBLTRANS <- bitr(k,fromType='ENSEMBLTRANS',toType=c('SYMBOL','ENTREZID'),OrgDb=orgdb,drop = TRUE)
    prefixes <- c('flank_','promoter_flank_'); targets <- c('SYMBOL','ENTREZID')
    sufixes <- c('txIds','geneIds'); keycols <- c('ENSEMBLTRANS','ENSEMBL')
    for (prefix in prefixes) { for (target in targets) { for (num in 1:2) {
        oldcol <- paste0(prefix,sufixes[num]); tmpcol <- paste0(prefix,target,sufixes[num]); df <- get(paste0('df_',keycols[[num]]))
        mcols(anot@anno)[[tmpcol]] <- sapply(mcols(anot@anno)[[oldcol]],get_flank_col,df=df,keycol=keycols[num],targetcol=target)
    }
    newcol <- paste0(prefix,target); tmpcols <- paste0(newcol,sufixes); print(newcol)
    mcols(anot@anno)[[newcol]] <- sapply(
        paste(mcols(anot@anno)[[tmpcols[[1]]]],mcols(anot@anno)[[tmpcols[[2]]]],sep=';'),del_NA_dup,USE.NAMES=FALSE)
    mcols(anot@anno)[[tmpcols[[1]]]] <- NULL; mcols(anot@anno)[[tmpcols[[2]]]] <- NULL
    }}
    print(paste0('########## 导出所有背景peak注释结果 ',name,'.annotation.xls.gz'))
    anot%>%as_tibble()%>%write_tsv(str_c(name,'.annotation.xls.gz'))
    return(anot)
}
plotAndMergeAnot <- function(name,peak,txdb,orgdb,anot.bg,flkdist=50000,promoter_upstream=2000,promoter_downstream=0){
    library(ChIPseeker)
    library(orgdb,character.only=TRUE)
    print('########## 对显著差异的peak进行注释')
    anot <- lapply(peak, annotatePeak,assignGenomicAnnotation=TRUE,
        tssRegion = c(-promoter_upstream, promoter_downstream),TxDb =txdb,annoDb=orgdb)
    tb.bg <- as_tibble(anot.bg)
    anot <- lapply(anot,function(i){
        tb.i <- as_tibble(i); tb <- left_join(tb.i[,1:3],tb.bg); selectcols <- !names(tb)%in%names(tb.i)
        mcols(i@anno) <- cbind(mcols(i@anno),tb[,selectcols])
        return(i)
    })
    anot$bg <- anot.bg
    print('########## 绘制peak分布图')
    pdf(str_c(name,'.coverage.distribution.pdf'),height=3,width=6)
        # print(covplot(peak))
        # print(plotAnnoPie(anot))
        print(plotAnnoBar(anot))
        # peak距离TSS的相对位置分布
        print(plotDistToTSS(anot))
    dev.off()
    # # 窗口结合图谱，需要资源较大
    # peakHeatmap(peak, TxDb = txdb, upstream = 1000, downstream = 1000)

    # 导出注释结果，这里的输出结果就不需要all了
    data.table::rbindlist(lapply(within(anot, rm(all,bg)),as.data.frame),fill=TRUE,idcol=F)%>%write_tsv(str_c(name,'.annotation.xls'))
    return(anot)
}
entrezidEnrich <- function(name,anot,orgdb,organism,flkdist=50000,promoter_upstream=2000,promoter_downstream=0){
    library(orgdb,character.only=TRUE)
    # 由于富集分析加入了pvalue和qvalue筛选，结果很可能是NULL，所以需要加入很多的异常处理。否则会报错停止
    # KEGG分析
    library(clusterProfiler)
    # gene选择nearest还是flank
    nearest<-lapply(within(anot, rm(all,bg)),function(i) as.data.frame(i)$ENTREZID)
    promoter.nearest<-lapply(within(anot, rm(all,bg)),function(i) as.data.frame(i)%>%
        dplyr::filter(distanceToTSS>=-promoter_upstream & distanceToTSS<=promoter_downstream)%>%.$ENTREZID)
    promoter.flank<-lapply(within(anot,rm(all,bg)),function(i)as.vector(unlist(strsplit(as.data.frame(i)$promoter_flank_ENTREZID,";"))))
    flank<-lapply(within(anot,rm(all,bg)),function(i)as.vector(unlist(strsplit(as.data.frame(i)$flank_ENTREZID,";"))))
    gene<-list(nearest,promoter.nearest,promoter.flank,flank)
    names(gene)<-c('nearest','promoter.nearest','promoter.flank',str_c('flank',flkdist))
    print('########## KEGG 富集分析开始')
    compKEGG <- lapply(gene,compareCluster,
        # geneCluster = gene,
        fun = "enrichKEGG",organism=organism,
        qvalueCutoff  = 0.05,pvalueCutoff  = 0.05,
        )
    # # 如果要自定义背景基因集，则使用universe参数，但是运算速度非常慢
    # bg.nearest<-as.data.frame(anot$bg)$ENTREZID
    # compKEGG$nearest_bg <- compareCluster(
    #     geneCluster = gene$nearest,
    #     fun = "enrichKEGG",organism=organism,
    #     qvalueCutoff  = 0.05,pvalueCutoff  = 0.05,
    #     universe = bg.nearest,
    #     )
    compKEGG <- compKEGG[!sapply(compKEGG, is.null)]  # 过滤掉没有结果的KEGG富集
    # 去除掉Description后面的宿主名称
    compKEGG<-lapply(compKEGG,function(x){
        if(!is.null(x)){ x@compareClusterResult$Description <- 
            str_replace(x@compareClusterResult$Description,'- Mus.*|- Homo.*','')
        };return(x)
    })
    saveRDS(compKEGG,file=str_c(name,'.KEGG.rds'))
    tb<-lapply(names(compKEGG),function(x)try(setReadable(compKEGG[[x]],OrgDb=orgdb,keyType="ENTREZID")%>%
        as_tibble()%>%write_tsv(str_c(name,'.KEGG.',x,'.xls'))))    # geneID转换成SYMBOL

    print('########## GO富集分析')
    # GO结果太多，简化方法：https://www.jianshu.com/p/cd30d9ad0129；https://github.com/ixxmu/mp_duty/issues/4773
    compGO <- lapply(gene,compareCluster,
        # geneCluster = gene,
        fun = "enrichGO",readable=T,OrgDb=orgdb,ont='ALL',
        qvalueCutoff  = 0.05,pvalueCutoff  = 0.05,
        )
    # # 如果要自定义背景基因集，则使用universe参数，但是运算速度非常慢
    # bg.nearest<-as.data.frame(anot$bg)$ENTREZID
    # compGO$nearest_bg <- compareCluster(
    #     geneCluster = gene$nearest,
    #     fun = "enrichGO",readable=T,OrgDb=orgdb,ont='ALL',
    #     qvalueCutoff  = 0.05,pvalueCutoff  = 0.05,
    #     universe = bg.nearest,
    #     )
    compGO <- compGO[!sapply(compGO, is.null)]  # 过滤掉没有结果的GO富集
    saveRDS(compGO,file=str_c(name,'.GO.rds'))
    compGO <- lapply(compGO,simplify,cutoff=0.7,by="p.adjust",select_fun=min)   # 简化GO结果，cutoff越大越精简
    tb<-lapply(names(compGO),function(x)try(as_tibble(compGO[[x]])%>%write_tsv(str_c(name,'.GO.',x,'.xls'))))
    # 绘图
    # library(enrichplot)
    # p <- pairwise_termsim(compKEGG)
    # treeplot(p)
    # cnetplot(compGO)
    pdf(str_c(name,'.KEGG.GO.pdf'))     # mapply这里面强制指定对象画图，会出现NULL Error是正常的
    tb<-mapply(function(x,y)try(print(dotplot(get(x)[[y]],showCategory=10,label_format=50,  # size='count',
        title=str_c(y,' ',str_sub(x,5)," Enrich")))),rep(c('compKEGG','compGO'),3),rep(names(gene),each=2))
    dev.off()
}

library(GenomicFeatures)
diff_peak.files <- str_split(diff_peak.files,'\\n')[[1]]
txdb <- makeTxDbFromGFF(gtf.file)
peak_counts_file <- str_c(dirname(dirname(diff_peak.files[[1]])),'/all/RAW_READS.peak_counts.xls.gz')
all.bg <- read_tsv(peak_counts_file)%>%filter(!str_detect(CHR,'_|V'))%>%makeGRangesFromDataFrame()
print('########## 开始注释背景peak')
anot.bg <- peakAnot(peak_counts_file,peak=all.bg,txdb,orgdb)
for (f in diff_peak.files){
    print(paste('######### 开始差异peak文件注释富集分析 ',basename(dirname(f))))
    name=str_c(dirname(f),'/peakAnnotationEnrich/',basename(f))
    if(!file.exists(str_c(name,'.annotation.xls'))){
        groupName <- paste0(strsplit(basename(dirname(f)), '-vs-')[[1]],'_enrich')
        groupName <- str_replace(groupName,'^LIB_FULL.|^RLE_BG.|^RLE_RIP.','')  # 要删掉前缀
        peak <- read_tsv(f)%>%filter(!str_detect(seqnames,'_|V'))%>%filter(abs(log2FoldChange)>1,FDR<0.05)%>%
            arrange(desc(log2FoldChange))%>%mutate(
                change=ifelse(log2FoldChange>0,groupName[1],groupName[2]))%>%rbind(.,mutate(.,change='all'))
        nrow_peak <- nrow(peak); if(nrow_peak<=0){next}     # peak数为0则跳过
        dir.create(dirname(name))
        # 创建GRanges对象，分成up、down和all三种
        peak <- peak%>%makeGRangesFromDataFrame(keep.extra.columns=T)%>%split(.$change)
        anot <- plotAndMergeAnot(name,peak,txdb,orgdb,anot.bg)
        if(nrow_peak<10 || nrow_peak>1e4 || enrich==''){next}      # peak10个以下或超过1万个，不做富集分析
        entrezidEnrich(name,anot,orgdb,organism)}
}
quit(save='no')



library(pathview)
data(gse16873.d)
pv.out <- pathview(gene.data = gse16873.d[, 1], pathway.id = "04110",
                   species = "hsa", out.suffix = "gse16873")
# 查看5个示例文件（Bioconductor自带示例数据）
files <- getSampleFiles()
basename(unlist(files))
# readPeakFile读取bed文件返回GRange对象,GRange是生信分析里约定俗成的一种结构
CBX6 <- readPeakFile(files[[4]])
CBX7 <- readPeakFile(files[[5]])
CBX7
# 直接画peak在各个染色体的分布图
covplot(CBX6)
# 多个文件同时作图（组成一个GRangeList对象），指定染色体及范围（MACS的BED第V5列为qvalue）
peaks <- list(CBX6 = CBX6, CBX7 = CBX7)
peak45 <- GenomicRanges::GRangesList(CBX6 = CBX6, CBX7 = CBX7)
covplot(peak45, weightCol = "V5", 
    chrs = c("chr17", "chr18"), xlim = c(4e7, 5e7)) +
    facet_grid(chr ~ .id)

# 画peak在某±1000的结合图谱和结合强度图
peakHeatmap(CBX6, TxDb = txdb, upstream = 1000, 
    downstream = 1000, color = rainbow(length(CBX6)))
# 也可以多个文件同时画
files=getSampleFiles()
peakHeatmap(files, TxDb=txdb, 
            upstream=3000, downstream=3000, 
            color=rainbow(length(files)))
# 以下是一步一步手动进行的方式：
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix <- getTagMatrix(f, windows=promoter)
tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="red")
# 画结合强度图：先直接使用bed文件
plotAvgProf2(CBX6, TxDb = txdb, upstream = 1000, 
    downstream = 1000, xlab = "Genomic Region (5'->3')", 
    ylab = "Read Count Frequency", conf = 0.95, resample = 1000)# 添加置信区间
# 手动进行：直接调用上面的tagMatrix
plotAvgProf(tagMatrix, xlim=c(-3000, 3000),
            xlab="Genomic Region (5'->3')", 
            ylab = "Read Count Frequency")
# peakHeatmap无法画多个，因此多个文件同时画需要采用分拆步骤：
promoter <- getPromoters(TxDb = txdb, upstream = 1000, downstream = 1000)
tagMatrixList <- lapply(peak45, getTagMatrix, windows = promoter)
tagHeatmap(tagMatrixList, xlim = c(-1000, 1000), color = "red")
# 添加置信区间
plotAvgProf(tagMatrixList, xlim = c(-1000, 1000), 
    conf = 0.95, resample = 500, facet = "row")

# https://bioconductor.org/packages/3.11/data/annotation/，官网上的物种不多
# 单个注释
CBX6anot <- annotatePeak(CBX6, tssRegion = c(-1000, 1000), 
    TxDb = txdb, annoDb = "org.Hs.eg.db")
# 多个同时注释（lapply对list的每个对象都执行同一个函数并返回结果list），verbose=FALSE不打印多余信息
peak45anot <- lapply(peak45, annotatePeak, tssRegion = c(-1000, 1000), 
    TxDb = txdb, annoDb = "org.Hs.eg.db", verbose = FALSE)
# 查看注释内容前三行
head(as.GRanges(CBX6anot))
head(lapply(peak45anot, as.GRanges))
# 转换为datafram形式保存
tmp <- as.data.frame(CBX6anot)
write.csv(tmp, "humanCBX6PeakAnot.csv")

# 非模式生物：通过makeTxDbFromGFF自己构建TxDb，也可以makeTxDbFromBiomart通过ensembl在线制作，还有其它函数
# https://bioconductor.org/packages/release/bioc/vignettes/GenomicFeatures/inst/doc/GenomicFeatures.pdf
# 注意：由gff构建而来的，后面注释的geneId会变成locus_tag，transcriptId会变成geneName，要做富集分析的话，需要使用自己根据locus_tag构建的GO和KEGG注释库
require(GenomicFeatures)
txdb <- makeTxDbFromGFF("~/project/refSeq/UCSC/mm39.ensGene.gtf")
columns(txdb)
bed <- readPeakFile("peaks.bed")
bed
anot <- annotatePeak(bed, tssRegion = c(-1000, 1000), TxDb = txdb)
as.data.frame(anot)%>%write_tsv("peakanot.xls")

# 画peak annotation构成，Pie图只能画单个，柱状图和TSS相对位置可以多个一起展示
vennpie(CBX6anot)
plotAnnoPie(CBX6anot)
library(ggimage)
upsetplot(CBX6anot, vennpie=TRUE)
plotAnnoBar(peak45anot)
# peak距离TSS的相对位置分布
plotDistToTSS(peak45anot, title = "Distribution of TF-binding loci relative to TSS")

# 不同bed文件基因交集韦恩图
# 先得到基因列表
genes <- lapply(peak45anot, function(i) 
    as.data.frame(i)$geneId)
names(genes)
vennplot(genes)

# 基因注释+富集分析
require(clusterProfiler)
files=getSampleFiles()
# 将bed文件读入（readPeakFile是利用read.delim读取，然后转为GRanges对象）
seq <- lapply(files, readPeakFile)

genes <- lapply(seq, function(i) 
    seq2gene(i, c(-1000, 3000), 3000, TxDb=txdb))
# 下面这一步直接运行会出错，https://cloud.tencent.com/developer/article/2032030
# 直接按教程也不对，修改以后如下
install.packages('R.utils')
R.utils::setOption("clusterProfiler.download.method",'wget')
cc <- compareCluster(geneClusters = genes, 
                    fun = "enrichKEGG", organism = "hsa")
dotplot(cc, showCategory = 10)
