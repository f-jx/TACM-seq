# https://www.cruk.cam.ac.uk/wp-content/uploads/2021/02/Quantitative-ChIPseq-Workshop-1.html
# https://github.com/hbctraining/In-depth-NGS-Data-Analysis-Course/blob/master/sessionV/lessons/08_diffbind_differential_peaks.md
# https://www.jianshu.com/p/b74c8077d893

# setwd('/thinker/glusterfs/home/QMLab-jxFang/project/tcga/all')
args <- c('diffbind_filter.csv','result/atac_filter/peak_differentialAnalysis/diffbind_filter.csv',16)
args <- commandArgs(trailingOnly = TRUE)
diffbind_input_file <- args[[1]]
contrast.group.file <- paste0('contrast.',diffbind_input_file)
result.dir <- args[[2]]
threads <- args[[3]]

library(DiffBind)
library(tidyverse)
library(pheatmap)

padj=0.05;fc=1
samples <- read.csv(diffbind_input_file,header=TRUE)    # %>%arrange(desc(Condition))
if(file.exists(contrast.group.file))contrastdf <- read.csv(contrast.group.file)
dir.create(result.dir,recursive=TRUE); setwd(result.dir); dir.create('all')
    
if (file.exists(str_c('dba_counts.RData'))){
    print(paste0('########## ','dba_counts.RData','文件存在，将直接加载'))
    counts <- dba.load(file='counts', dir='.', pre='dba_', ext='RData')
} else {
    print('########## dba读取samplesheet文件')
    peaks <- dba(sampleSheet=samples)
    peaks$config$RunParallel <- TRUE
    peaks$config$cores <- threads
    # paste0(names(peaks),collapse=', ')  # merged存储merge的peaks，但seqname变成了编号，config是相关的设置参数
    # peaks$peaks[[1]]%>%as_tibble()  # "Chr, Start, End, Score"

    # # generate a greylist，添加peak黑名单
    # data(tamoxifen_greylist)
    # names(tamoxifen.greylist)
    # names(tamoxifen.greylist$controls)
    # tamoxifen.greylist$master
    # names(tamoxifen.greylist$controls)
    # # apply the greylist
    # peaks <- dba.blacklist(peaks, greylist=tamoxifen.greylist)

    # 计算这些每一个peak上不同样品的重叠reads数量，无论该样品中是否有该peak，生成矩阵
        # 参数：peaks限定统计区间，summits会按照峰顶重新定义peak区间（默认=200,以峰顶为中心上下游共401bp），
        # filter过滤至少1个样品RPKM超过1，bRemoveDuplicates默认不删除dup reads，
        # mapQCth忽略低质量比对，默认15，bUseSummarizeOverlaps使用比默认更标准的计数
        # bScaleControl和bSubControl和ChIP seq有关
    print('########## 执行dba.count计数')
    counts <- dba.count(peaks,summits=FALSE,filter=1,bRemoveDuplicates=FALSE,     # ,summits=200
        bScaleControl=FALSE,mapQCth=30,bUseSummarizeOverlaps=TRUE,bParallel=TRUE)
    # counts中多了两列，Reads显示样品总reads数，及Fraction of Reads in Peaks，FRiP
    # paste0(names(counts),collapse=', ') # 比peaks多出了score，SN，maxFilter，filterFun，minCount, summits
    # 这时的merged已经是调整后的peaks，peaks中多出了4列，RPKM,Reads,cRPKM,cREADS（c可能代表ChIPseq）
    # 所有的样品都针对minOverlap=2的peak计算了Reads数，行数一样，仔细计算一下，发现RPKM有除以peak宽度
    # counts$peaks[[1]]%>%as_tibble() # "Chr, Start, End, Score, RPKM, Reads, cRPKM, cReads"; peak调整过后的结果
    dba.save(counts,file='counts')
    # 保存数据表格
    dba.show(peaks)%>%  # 保存样品分组信息
        left_join(dba.show(counts),by=c('ID','Condition'),suffix=c('_peak','_count'))%>%select(!contains('Replicate'))%>%
        left_join(mutate(peaks$samples,ID=SampleID,SampleID=NULL,Replicate=NULL))%>%write_tsv('all/samples.xls')
    # # consensus peak set，只是最初的merge peak结果(minOverlap=2)，之后会调整(recenter 401bp，filter等)
    # dba.peakset(peaks, bRetrieve=TRUE)[,0]%>%as_tibble()%>%transmute(seqnames,start=start-1,end)%>%
    #     write_tsv("all/consensus.peak.minOverlap2.bed",col_names=FALSE)
    # 绘图展示结果
    pdf('all/overlap.rate.pdf')
        olap.rate <- dba.overlap(peaks, mode=DBA_OLAP_RATE)
        plot(olap.rate, xlab="Overlapping samples", ylab="Overlapping peaks", type="b")
        text(olap.rate, as.character(olap.rate),pos=3,col='red')
    dev.off()
}

print('########## 执行dba.normalize标准化')
# 数据归一化，比counts多了一列norm，记录归一化方法。这里选择两种最常用的，具体方案查看Diffbind的Normalize一章最后结论。
# lib + full是按测序文库深度做标准化不做任何假设得到最真实的样品结果；
# rle + bg假设大多数feature是没有显著改变的，用除了peak以外的背景区域来计算size factor。
counts_norm.LIB_FULL <- dba.normalize(counts,method=DBA_DESEQ2,normalize=DBA_NORM_LIB,library=DBA_LIBSIZE_FULL,background=FALSE)
if (file.exists(str_c('dba_RLE_BG.RData'))){
    counts_norm.RLE_BG <- dba.load(file='RLE_BG', dir='.', pre='dba_', ext='RData')
} else {
    counts_norm.RLE_BG <- dba.normalize(counts,method=DBA_DESEQ2,normalize=DBA_NORM_RLE,library=DBA_LIBSIZE_BACKGROUND,background=TRUE)    # BiocManager::install("csaw",ask=FALSE) background需要安装csaw
    dba.save(counts_norm.RLE_BG,file='RLE_BG')
}
counts_norm.RLE_RIP <- dba.normalize(counts,method=DBA_DESEQ2,normalize=DBA_NORM_RLE,library=DBA_LIBSIZE_PEAKREADS,background=FALSE)
# paste0(names(counts_norm),collapse=', ')    # 多了一个norm值

print('########## 输出样本的LibSize和FRiP')
info<-dba.show(counts_norm.LIB_FULL); norm <- dba.normalize(counts_norm.LIB_FULL, bRetrieve=TRUE)
normlibs<-cbind(sample=info$ID,FullLibSize=norm$lib.sizes,FRiP=info$FRiP,PeakReads=round(info$Reads * info$FRiP))
normlibs%>%as_tibble()%>%write_tsv('all/FullLibSize.FRiP.xls')

print('########## 保存rawCounts、不同标准化后的矩阵，绘制相关性热图和PCA图')   # dba.count(counts,peaks=NULL,score=DBA_SCORE_READS) # 获取原始Counts表
tb <- dba.peakset(dba.count(counts,peaks=NULL,score=DBA_SCORE_READS),bRetrieve=TRUE,DataType=DBA_DATA_FRAME)
colnames(tb)[-1:-3] <- counts$samples$SampleID  # 上一步dba.peakset可能会自动修改样品列名，这里要改回来
tb%>%write_tsv('all/RAW_READS.peak_counts.xls.gz'); rm(tb)
# 对所有样本绘制相关性热图，然后选择其中的一列条件，画出PCA聚类图，这里选condition，也可以是DBA_TISSUE。
pdf(paste0("all/RAW_READS.correlation.pca.pdf"))
    dh <- dba.plotHeatmap(counts,score=DBA_SCORE_READS,attributes=DBA_CONDITION,ColAttributes=DBA_CONDITION)
    # dh %>% write_tsv(paste0('all/',norm,'.correlation.xls'))    # dh是个double matrix，不能直接保存为xls
    dba.plotPCA(counts,score=DBA_SCORE_READS,attributes=DBA_CONDITION,vColors=ggsci::pal_d3("category20")(20))
dev.off()
tmp <- mapply(function(dba_obj,norm){
    # 用dba.peakset(counts_norm,bRetrieve=TRUE)提取，其实就是peaks$peaks[[1]]%>%as_tibble()中的score列
    tb <- dba.peakset(dba_obj,bRetrieve=TRUE,DataType=DBA_DATA_FRAME)
    colnames(tb)[-1:-3] <- dba_obj$samples$SampleID  # 上一步dba.peakset可能会自动修改样品列名，这里要改回来
    tb%>%write_tsv(paste0('all/',norm,'.peak_counts.xls.gz'))
    pdf(paste0("all/",norm,".correlation.pca.pdf"))
        dh <- dba.plotHeatmap(dba_obj,score=DBA_SCORE_NORMALIZED,attributes=DBA_CONDITION,ColAttributes=DBA_CONDITION)
        # dh %>% write_tsv(paste0('all/',norm,'.correlation.xls'))    # dh是个double matrix，不能直接保存为xls
        dba.plotPCA(dba_obj,score=DBA_SCORE_NORMALIZED,attributes=DBA_CONDITION,vColors=ggsci::pal_d3("category20")(20))
    dev.off()
    },list(counts_norm.LIB_FULL,counts_norm.RLE_BG,counts_norm.RLE_RIP),c('LIB_Full','RLE_BG','RLE_RIP'))
rm(tmp)
hmap<-colorRampPalette(c('#2878B5','#C82423'))(n=200)
# heatmap默认只画1000行，且无法取消聚类，暂时不用
# pdf('all/heatmap.pdf')
# p.all<-dba.plotHeatmap(counts,correlations=FALSE,score=DBA_SCORE_RPKM,scale='row',colScheme=hmap)
# dev.off()

# 判断是否有差异分析分组文件存在
if(!exists('contrastdf')){quit(save='no')}
print('########## contrastdf文件存在，将执行差异分析')

dba_analyze_plot <- function(counts_norm,norm,contrast.df=contrastdf){
print('########## 执行dba.contrast、dba.analyze创建对比组并进行差异分析')
# 先建一个不带contrast的model，然后逐一添加contrast，若不指定，会默认对所有分组进行两两比较
model <- dba.contrast(counts_norm,design='~Condition',minMembers = 2)
for(i in 1:nrow(contrast.df)){
    treatment=contrast.df[i,1];ctrl=contrast.df[i,2]
    # if (nrow(counts$samples[counts$samples$Condition==treatment,])<2 | nrow(counts$samples[counts$samples$Condition==ctrl,])<2){next}   # 样品不足跳过
    contrast <- c('Condition',treatment,ctrl)
    print(contrast)
    model<-dba.contrast(model,contrast=contrast)    # 添加要对比的组别
}    
# paste0(names(model),collapse=', ')              # model中多出，meta，design，contrasts等
# model$contrast                                  # contrast中有许多内容，可以看出设置contrast的方式
# 添加完contrast后，统一进行分析，再按类别分开导出
result <- dba.analyze(model)
# paste0(names(result),collapse=', ')

# 合并所有上下调peaks
merge.up<-data.frame(matrix(nrow=0,ncol=3));colnames(merge.up)<-c('seqnames','start','end'); merge.dw<-merge.up
# 准备heatmap调色板
color<-c(colorRampPalette(c('#2878B5','white'))(n=401),colorRampPalette(c('white','#C82423'))(n=400))
bk<-c(seq(-2,0,by=0.005),seq(0.005,2,by=0.005))
tb <- dba.peakset(result,bRetrieve=TRUE,DataType=DBA_DATA_FRAME);colnames(tb)[-1:-3] <- result$samples$SampleID   # 后面绘制heatmap要用到，但样品列名要改回来
print('########## 查看需要执行差异分析的组别')
contrasts<-dba.show(result,bContrasts=TRUE)
print(contrasts)
for (i in 1:nrow(contrasts)){
    # i=5
    factor = contrasts[i,1]; treat = contrasts[i,2]; ctrl = contrasts[i,4]
    print(paste0('########## 生成差异分析结果 ',norm,'.',treat,'-vs-',ctrl))
    dir <- paste0(norm,'.',treat,'-vs-',ctrl); dir.create(dir)
    # 生成DESeq2分析后的结果报告，其中每一列的意义查看网站描述；这里有个非常严重的问题容易混淆，如果在参数中设置了fold值，那么所有的FDR和P值都会在fold过滤的基础上重新计算，和不设置fold的完全不一样！
    # rep<-dba.report(result, contrast= i,th=padj,fold=fc)%>%as_tibble();rep%>%arrange(p.value)%>%rename(log2FoldChange=Fold)%>%write_tsv(str_c(dir,"/diff_peaks.th0.05.xls"))
    rep<-dba.report(result, contrast= i, th = 1)%>%as_tibble()    # 参数th 默认选择FDR<0.05的，所以可能为空；先把条件放宽，之后再筛选

    print('########## draw PCA MA and Volcano')
    pdf(str_c(dir,'/PCA.Volcano.pdf'))
        my_cols <- c("#C82423", "gray", "#2878B5"); names(my_cols) <- c("up", "none", "down")
        data <- rep%>%mutate(significant='none')
        try(dba.plotPCA(result,contrast=i,th=padj,attributes=DBA_CONDITION,label=NULL,vColors=ggsci::pal_d3("category20")(20)))
        try(dba.plotMA(result,contrast=i,th=padj,bUsePval=FALSE,dotSize=0.5,bSignificant=FALSE))    # ,fold=fc 注意，加了fold后，p值和FDR会重新计算！
        # dba.plotVolcano(result,contrast=i,th=padj,fold=fc)
        # 手动绘制MA图和火山图，可以考虑用这个包ggpubr::ggmaplot
        p <- ggplot(data,aes(Conc, Fold)) + geom_point(aes(color=significant),alpha=0.6)+
            ylim(-max(abs(data$Fold)),max(abs(data$Fold)))+
            labs(title='MA plot', x='log concentration', y=expression(log[2](FC)))+
            scale_color_manual(values=my_cols)+     # ,drop=FALSE 无论有没有这个level，都会在图例中显示出来
            geom_hline(yintercept=c(-1,1),linetype=4)+
            theme_bw()+theme(panel.grid=element_line(color='gray97'),
                text=element_text(size=16,color='black'),
                axis.text=element_text(color='black'),
                plot.title=element_text(hjust = 0.5,size=18),
                axis.text.x=element_text(angle=90,hjust=1))
        print(p)
        data[data$Fold < -fc & data$FDR < padj, c("significant")] = "down"
        data[data$Fold > fc & data$FDR < padj, c("significant")] = "up"
        # data$significant <- factor(data$significant, levels = c("down", "none", "up"))
        p <- ggplot(data,aes(Fold, -1*log10(FDR))) + geom_point(aes(color=significant))+
            xlim(-max(abs(data$Fold)),max(abs(data$Fold)))+ylim(0,max(-log10(data$FDR)))+
            labs(title='Volcanoplot', x=expression(log[2](FC)), y=expression(-log[10](p.adjust)))+
            scale_color_manual(values=my_cols)+     # ,drop=FALSE 无论有没有这个level，都会在图例中显示出来
            geom_hline(yintercept=-log10(padj),linetype=4)+
            geom_vline(xintercept=c(-1,1),linetype=4)+
            theme_bw()+theme(panel.grid=element_line(color='gray97'),
                text=element_text(size=16,color='black'),
                axis.text=element_text(color='black'),
                plot.title=element_text(hjust = 0.5,size=18),
                axis.text.x=element_text(angle=90,hjust=1))
        print(p)
        # 如果log2FC差异非常小，会导致所有的FDR接近1，画出来的火山图很奇怪，所以改用P值绘图
        data[data$Fold < -fc & data$p.value < padj, c("significant")] = "down"
        data[data$Fold > fc & data$p.value < padj, c("significant")] = "up"
        # data$significant <- factor(data$significant, levels = c("down", "none", "up"))
        p <- ggplot(data,aes(Fold, -1*log10(p.value))) + geom_point(aes(color=significant))+
            xlim(-max(abs(data$Fold)),max(abs(data$Fold)))+ylim(0,max(-log10(data$p.value)))+
            labs(title='Volcanoplot', x=expression(log[2](FC)), y=expression(-log[10](Pvalue)))+
            scale_color_manual(values=my_cols)+     # ,drop=FALSE 无论有没有这个level，都会在图例中显示出来
            geom_hline(yintercept=-log10(padj),linetype=4)+
            geom_vline(xintercept=c(-1,1),linetype=4)+
            theme_bw()+theme(panel.grid=element_line(color='gray97'),
                text=element_text(size=16,color='black'),
                axis.text=element_text(color='black'),
                plot.title=element_text(hjust = 0.5,size=18),
                axis.text.x=element_text(angle=90,hjust=1))
        print(p)
    dev.off()
    rep%>%arrange(p.value)%>%rename(log2FoldChange=Fold)%>%write_tsv(str_c(dir,"/diff_peaks.xls"))
    # 保存bed文件，用于后续作图，注意保存成bed文件需要把开始位置-1，而xls不需要
    up<-rep%>%filter(Fold>fc,FDR<padj)%>%transmute(seqnames,start,end)
    dw<-rep%>%filter(Fold< -fc,FDR<padj)%>%transmute(seqnames,start,end)
    if(nrow(up)>0)up%>%transmute(seqnames,start=start-1,end)%>%write_tsv(str_c(dir,'/up.bed'),col_names=FALSE)
    if(nrow(dw)>0)dw%>%transmute(seqnames,start=start-1,end)%>%write_tsv(str_c(dir,'/down.bed'),col_names=FALSE)
    share<-rep%>%filter(abs(Fold)>fc,FDR<padj)%>%transmute(seqnames,start=start-1,end)
    if(nrow(share)>0)share%>%write_tsv(str_c(dir,'/share.bed'),col_names=FALSE)
    merge.up<-rbind(merge.up,up)    # 将上下调peaks合并入merge表格
    merge.dw<-rbind(merge.dw,dw)
    if(nrow(share)<5){print('########## 两组差异paak<5，跳过'); file.create(str_c(dir,'/差异peak数小于5.tmp')); next}
    print('########## 生成热图表格并绘制差异peaks热图')
    hmp<-rep%>%filter(abs(Fold)>fc,FDR<padj)%>%arrange(desc(Fold))%>%
        transmute(CHR=seqnames,START=start,END=end)%>%left_join(tb)%>%mutate(peak=str_c(CHR,'_',START,'_',END))%>%
        data.frame(row.names=.$peak,check.names=F)%>%mutate(peak=NULL,CHR=NULL,START=NULL,END=NULL)
    anot_col<-samples%>%as_tibble()%>%transmute(SampleID,Condition)%>%  # arrange(Condition)%>%
        data.frame(row.names=.$SampleID,check.names=F)%>%mutate(SampleID=NULL)
    ctrls <- result$samples%>%filter(Condition==ctrl)%>%transmute(SampleID)%>%.[[1]]
    treats <- result$samples%>%filter(Condition==treat)%>%transmute(SampleID)%>%.[[1]]
    # 画单独的contrast图
    pheatmap(hmp[,c(ctrls,treats)],cluster_cols=F,cluster_rows=F,angle_col="90",scale='row',
        breaks=bk,color=color,annotation_col=anot_col[c(ctrls,treats),,drop=F],gaps_row=nrow(up),
        gaps_col=cumsum(rle(anot_col[c(ctrls,treats),,drop=F]$Condition)$lengths),show_rownames = F,
        main=str_c("Significantly DEPs expression(row zscore)",nrow(up),nrow(dw),sep=' '),width=7,height=7,
        filename=str_c(dir,'/DEPs.',nrow(up),'_',nrow(dw),'.heatmap.pdf'))
    hmp[,c(ctrls,treats)]%>%as_tibble(rownames='peak')%>%
        write_csv(str_c(dir,'/DEPs.',nrow(up),'_',nrow(dw),'.heatmap.csv'))

    # pdf(str_c(dir,'/heatmap.pdf'))
    # p.vs<-dba.plotHeatmap(result,contrast=i,correlations=FALSE,scale='row',colScheme=hmap)
    # dev.off()
    # 画deeptools类似的热图，但不会调整，暂时不用
    # require(profileplyr)
    # repList<-GRangesList(
    #     Gain=rep[rep$Fold>1,],
    #     Share=rep[abs(rep$Fold)<=1,],
    #     Loss=rep[rep$Fold< -1,]
    #     )
    # sampList<-list(result$masks[treat][[1]],result$masks[ctrl][[1]])
    # names(sampList)<-c(treat,ctrl)
    # b<-dba.plotProfile(result,contrasts=1,sites=repList,
    #     samples=sampList,
    #     merge=NULL,
    #     # labels=list(samples=c(treat,ctrl)),
    #     # bin_size=200,distanceAround=5000,
    #     )
    # pdf(str_c(dir,'.plotProfile.pdf'))
    # dba.plotProfile(b)
    # dev.off()
}
print('########## 所有contrast merge起来的图')
distinct(merge.up)%>%transmute(seqnames,start=start-1,end)%>%write_tsv(paste0('all/',norm,'.merge.up.bed'),col_names=FALSE)
distinct(merge.dw)%>%transmute(seqnames,start=start-1,end)%>%write_tsv(paste0('all/',norm,'.merge.down.bed'),col_names=FALSE)
merge<-distinct(rbind(distinct(merge.up),distinct(merge.dw)))
if(nrow(merge)<5){return()}
hmp<-merge%>%transmute(CHR=seqnames,START=start,END=end)%>%left_join(tb)%>%
    mutate(peak=str_c(CHR,'_',START,'_',END))%>%
    data.frame(row.names=.$peak,check.names=F)%>%mutate(peak=NULL,CHR=NULL,START=NULL,END=NULL)
anot_col<-samples%>%as_tibble()%>%transmute(SampleID,Condition)%>%  # arrange(Condition)%>%
    data.frame(row.names=.$SampleID,check.names=F)%>%mutate(SampleID=NULL)
pheatmap(hmp[,rownames(anot_col)],cluster_cols=F,cluster_rows=F,angle_col="90",scale='row',breaks=bk,color=color,
    annotation_col=anot_col,gaps_col=cumsum(rle(anot_col$Condition)$lengths),gaps_row=nrow(distinct(merge.up)),
    main=str_c("All significantly DEPs(zscore by row)",nrow(merge.up),nrow(merge.dw),sep=' '),width=10,height=10,show_rownames = F,
    filename=paste0('all/',norm,'.merge.heatmap.pdf'))
hmp%>%as_tibble(rownames='peak')%>%write_csv(paste0('all/',norm,'.merge.heatmap.csv'))
}
dba_analyze_plot(counts_norm.LIB_FULL,'LIB_FULL')
dba_analyze_plot(counts_norm.RLE_BG,'RLE_BG')
dba_analyze_plot(counts_norm.RLE_RIP,'RLE_RIP')
quit(save='no')



# 相关的绘图：
dba.plotMA(model)
dba.plotVolcano(model)
plot(model, contrast=1)
dba.plotPCA(model, contrast=1, label=DBA_TISSUE)
hmap <- colorRampPalette(c("red", "black", "green"))(n = 13)
readscores <- dba.plotHeatmap(model, contrast=1, correlations=FALSE,
    scale="row", colScheme = hmap)

# 多因素设计
model <- dba.contrast(model, design="~Tissue + Condition")
dba.show(model,bDesign=TRUE)
model <- dba.analyze(model)
dba.plotMA(model, contrast=1)

# 归一化方法
dba.plotMA(model,bNormalized=FALSE, th=0, sub="Non-Normalized")
dba.plotMA(model,bNormalized=TRUE, sub="Normalized: Library Size")
model <- dba.normalize(model, normalize="RLE", background=TRUE)
model <- dba.analyze(model)
dba.plotMA(model,bNormalized=TRUE, sub="Normalized: RLE [Background]")
