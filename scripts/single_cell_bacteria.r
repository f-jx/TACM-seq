library(tidyverse)
library(stringr)
# setwd('~/project/20230925_scATAC')
args<-c('sample.rename.group.csv',
    'merge/HSHF.readsCount_barcode.gz',            # 获得细胞的totalReads
    'kraken2/hg38_2Bc_StmFuso/HSHF.k2out.readID_txID.gz',      # 获得细胞kraken2分类结果
    'merge/kraken2.hg38_2Bc_StmFuso.taxon.xls',                       # 根据txID获得分类层级
    'bowtie2/host/HSHF/outs/filtered_tf_bc_matrix/barcodes.tsv.gz',    # 获得被保留下的细胞
    'singlecell/microbes',
    'Species','Salmonella enterica,Fusobacterium nucleatum')                          # 需要统计的层级和对应名称
args <- commandArgs(trailingOnly=TRUE)
sample.group <- read_csv(args[1])
sample <- sample.group%>%filter(sample==basename(dirname(dirname(dirname(args[5])))))%>%.$newname;sample
k2db <- str_replace_all(args[4],'.*kraken2.|.taxon.xls','');k2db

dir.create(args[6],recursive=TRUE)
if (!file.exists(str_c(args[6],'/',sample,'.',k2db,'.xls.gz'))){
    # 从cellranger bam文件统计所有barcode的reads数
    scbam<-read_tsv(args[2],col_names=c('totalReads','corrected_barcode'))
    # 统计unmap的每个readsID所属barcode，注意：有些unmap readsID在进入kraken2之前就被fastp过滤掉，因此没有对应txID信息
    unmap<-read_tsv(str_replace(args[2],'readsCount_barcode','unmapReadID_barcode'),
        col_names=c('readsID','corrected_barcode'))%>%mutate(readsID=str_replace_all(readsID,'/1$|/2$|@1$|@2$',''))%>%  # 有些空间组学的readsID后面会带双端的标识，需要去除
        group_by(corrected_barcode)%>%mutate(unmapReads=n())
    taxon<-read_tsv(args[4])%>%mutate(DBMinimizers=NULL,DBDistinctMinimizers=NULL)
    # 读取keep的barcode
    keep_cell<-read_tsv(args[5],col_names=c('corrected_barcode'))%>%mutate(group='keep')
    # 读取unmap每个readsID的txID，连接所属barcode并统计每个barcode，每个txID分配的kraken2 reads数
    k2out<-read_tsv(args[3],col_names=c('readsID','txID'))%>%
        full_join(unmap)%>%mutate(readsID=NULL)%>%
        group_by(corrected_barcode,txID,unmapReads)%>%
        summarise(kraken2Reads=n())%>%ungroup()     # txID为NA的kraken2Reads是被fastp filter掉的
    # 连接上述表格，barcode，totalReads，mappedReads，txID Reads(未作层级累加)，taxon信息
    k2out<-k2out%>%right_join(scbam,by='corrected_barcode')%>%
        mutate(mappedReads=ifelse(is.na(unmapReads),        # 计算mappedReads
            totalReads,totalReads-unmapReads),unmapReads=NULL)%>%
        left_join(keep_cell,by='corrected_barcode')%>%mutate(group=ifelse(is.na(group),
            ifelse(is.na(corrected_barcode),'no_tag','filtered'),group))%>%     # 添加细胞标签
        left_join(taxon,by='txID')      # 添加kraken2分类信息
    # 保存原始统计结果
    k2out%>%arrange(group,corrected_barcode)%>%
        write_tsv(str_c(args[6],'/',sample,'.',k2db,'.xls.gz'))
} else { k2out<-read_tsv(str_c(args[6],'/',sample,'.',k2db,'.xls.gz')) }



# 计算microbereads
k2out <- k2out%>%group_by(corrected_barcode,totalReads)%>%
    filter(Domain%in%c("Bacteria","Archaea","Viruses") | Kindom=="Fungi")%>%
    summarise(MicrobeReads=sum(kraken2Reads,na.rm=TRUE))%>%
    mutate(MicrobeReads=ifelse(is.na(MicrobeReads),0,MicrobeReads),
        lg2a1_microbeReads=log2(1+MicrobeReads),
        MicrobeRPM=1e6*MicrobeReads/totalReads,
        lg2a1_microbeRPM=log2(1+1e6*MicrobeReads/totalReads))%>%right_join(k2out)
# 统计不同group的reads数
group.sum<-k2out[,c('corrected_barcode','group','totalReads','mappedReads','MicrobeReads')]%>%
    distinct()%>%group_by(group)%>%summarise(
        Barcodes_counts=n(),TotalReads=sum(totalReads,na.rm=TRUE),
        MappedReads=sum(mappedReads,na.rm=TRUE),
        MicrobeReads=sum(MicrobeReads,na.rm=TRUE))

# 计算特别species的reads数量
summarise_reads <- function(S){
    k2out <<- k2out%>%filter(.data[[args[7]]]==S)%>%group_by(corrected_barcode,totalReads)%>%
        summarise(Reads=as.double(sum(kraken2Reads,na.rm=TRUE)))%>%
        mutate(Reads=ifelse(is.na(Reads),0,Reads),lg2a1_Reads=log2(1+Reads),
        RPM=1e6*Reads/totalReads,lg2a1_RPM=log2(1+1e6*Reads/totalReads))%>%right_join(k2out)    # 这里右连接，会引入新的NA值，需要修改
    group.sum <<- left_join(group.sum,k2out[,c('corrected_barcode','group','Reads')]%>%
        distinct()%>%group_by(group)%>%summarise(Reads=sum(Reads,na.rm=TRUE)))
    names(k2out)[names(k2out)=='Reads'] <<- paste0(substr(args[7],1,1),'_',gsub('[a-z]* ','.',S),'Reads')
    names(k2out)[names(k2out)=='RPM'] <<- paste0(substr(args[7],1,1),'_',gsub('[a-z]* ','.',S),'RPM')
    names(k2out)[names(k2out)=='lg2a1_Reads'] <<- paste0('lg2a1_',substr(args[7],1,1),'_',gsub('[a-z]* ','.',S),'Reads')
    names(k2out)[names(k2out)=='lg2a1_RPM'] <<- paste0('lg2a1_',substr(args[7],1,1),'_',gsub('[a-z]* ','.',S),'RPM')
    names(group.sum)[names(group.sum)=='Reads'] <<- paste0(substr(args[7],1,1),'_',gsub('[a-z]* ','.',S),'Reads')
}

scatter_plot <- function(var){
    library(cowplot)
    library(ggpubr)
    library(scales)
    # Main plot
    dt<-k2out[,c('corrected_barcode','totalReads','mappedReads','group',var)]%>%
        filter(group!='no_tag')%>%distinct()%>%mutate(mappedReads=as.double(mappedReads),
        alpha=ifelse(group=='filtered',0.3,0.7))
    dt[[var]]<-as.double(dt[[var]])

    # reads，dt[[var]]返回向量，dt[var]返回列表，0.5大概是10的-0.3
    dt[dt$mappedReads==0,'mappedReads']<-0.5
    dt[is.na(dt[[var]]),var]<-0.5
    value=log10(max(dt$mappedReads))+0.01
    value=6
    dt<-dt%>%mutate(x=log10(mappedReads),y=log10(.data[[var]]))
    pmain<-dt%>%ggplot(aes(x=x,y=y,color=group))+
        geom_point(alpha=0.3) + color_palette("jco")+   # 这里取-0.45~value之间，保证囊括所有点
        scale_x_continuous(limits=c(-0.45,value),breaks=c(0:floor(value)),labels=math_format(10^.x))+
        scale_y_continuous(limits=c(-0.45,value),breaks=c(0:floor(value)),labels=math_format(10^.x))+
        labs(x='cellrangerMapped reads counts',y=str_c(var,' reads counts in cellrangerUnmap'))+
        annotation_logticks(sides="trbl") + theme(text=element_text(size=14))
    xdens <- axis_canvas(pmain, axis="x") +
        geom_density(data=dt, aes(x=x, fill=group),
            position='stack',bw=0.1,color=NA,alpha=0.6, linewidth=0.2) +    # color=NA不显示轮廓线
        fill_palette("jco")
    ydens <- axis_canvas(pmain, axis="y", coord_flip=TRUE) +
        geom_density(data=dt, aes(x=y, fill=group),
            position='stack',bw=0.1,
            color=NA,alpha=0.6, linewidth=0.2) +
        coord_flip() + fill_palette("jco")
    # xdens <- axis_canvas(pmain, axis="x") +
    #     geom_histogram(data=dt, aes(x=x, fill=group),color='black',
    #         breaks=c(-0.5,-0.23,0,log10(2),log10(3),log10(5),log10(10),log10(20),log10(50),log10(100),log10(1e3),log10(1e6))) +
    #     scale_y_continuous(trans=log1p_trans())+fill_palette("jco")
    # ydens <- axis_canvas(pmain, axis="y", coord_flip=TRUE) +
    #     geom_histogram(data=dt, aes(x=y, fill=group),color='black',
    #         breaks=c(-0.5,-0.23,0,log10(2),log10(3),log10(5),log10(10),log10(20),log10(50),log10(100),log10(1e3),log10(1e6))) +
    #     scale_y_continuous(trans=log1p_trans())+
    #     coord_flip() + fill_palette("jco")
    p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position="top")
    p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position="right")
    pdf(str_c(args[6],'/',sample,'.',k2db,'.',var,'.reads.pdf'),width=8,height=7)
    print(ggdraw(p2))
    dev.off()

    # percent
    dt[as.vector(dt[var]==0.5),var]<-0
    dt<-dt%>%mutate(x=log10(mappedReads),y=.data[[var]]/totalReads)
    pmain<-dt%>%ggplot(aes(x=x,y=y,color=group,alpha=0.5))+
        geom_point(alpha=0.3) + color_palette("jco")+
        scale_x_continuous(limits=c(-0.45,value),breaks=c(0:floor(value)),labels=math_format(10^.x))+
        scale_y_continuous(limits=c(-0.01,1.01),labels=percent_format(accuracy=1))+
        labs(x='cellrangerMapped reads counts',y=str_c(var,' / (cellrangerMaped + unmap) percent'))+
        annotation_logticks(sides="tb") +theme(text=element_text(size=14))
    # # Marginal densities along x axis
    xdens <- axis_canvas(pmain, axis="x") +
        geom_density(data=dt, aes(x=x, fill=group),
            color=NA,alpha=0.6, linewidth=0.2) + fill_palette("jco")
    # Marginal densities along y axis
    # Need to set coord_flip=TRUE, if you plan to use coord_flip()
    ydens <- axis_canvas(pmain, axis="y", coord_flip=TRUE)+
        geom_density(data=dt, aes(x=y, fill=group),
            color=NA,alpha=0.6, linewidth=0.2) +
        coord_flip() + fill_palette("jco")
    # xdens <- axis_canvas(pmain, axis="x") +
    #     geom_histogram(data=dt, aes(x=x, fill=group),color='black',
    #         breaks=c(-0.5,-0.23,0,log10(2),log10(3),log10(5),log10(10),log10(20),log10(50),log10(100),log10(1e3),log10(1e6))) +
    #     scale_y_continuous(trans=log1p_trans())+
    #     fill_palette("jco")
    # ydens <- axis_canvas(pmain, axis="y", coord_flip=TRUE) +
    #     geom_histogram(data=dt, aes(x=y, fill=group),color='black',
    #         breaks=c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1)) +
    #     scale_y_continuous(trans=log1p_trans())+
    #     coord_flip() + fill_palette("jco")
    p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position="top")
    p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position="right")
    pdf(str_c(args[6],'/',sample,'.',k2db,'.',var,'.percent.pdf'),width=8,height=7)
    print(ggdraw(p2))
    dev.off()
}
# 绘制分箱图
plot_bin <- function(vars){
    df <- k2out[,c('corrected_barcode','totalReads','mappedReads','group',vars)]%>%distinct()%>%
        filter(group!='no_tag')%>%pivot_longer(cols=vars,names_to='type',values_to='reads')%>%
        mutate(reads=ifelse(is.na(reads),0,reads),RPM=1e6*reads/totalReads)%>%
        mutate(bins_reads=cut(reads,breaks=c(0,1,2,3,5,10,20,50,100,200,500,1e3,Inf),right=F),
            bins_RPM=cut(RPM,breaks=c(0,10,20,50,100,200,500,1e3,5e3,1e4,Inf),right=F))
    df_reads <- df%>%group_by(group,type,bins_reads)%>%summarise(n=n())%>%mutate(prop=n/sum(n))
    df_RPM <- df%>%group_by(group,type,bins_RPM)%>%summarise(n=n())%>%mutate(prop=n/sum(n))
    p_reads<-df_reads%>%ggplot(aes(y=bins_reads,x=n))+geom_col(position='dodge',color='gray',fill='lightblue')+
        geom_text(aes(label=paste0(n,'\n',round(prop,2))),position=position_dodge(0.9),size=2.7)+
        labs(title=sample,x='cell number',y='bacteria readsCount bins')+
        facet_grid(type~group,scale='free')+theme_classic()
    p_RPM<-df_RPM%>%ggplot(aes(y=bins_RPM,x=n))+geom_col(position='dodge',color='gray',fill='lightblue')+
        geom_text(aes(label=paste0(n,'\n',round(prop,2))),position=position_dodge(0.9),size=2.7)+
        labs(title=sample,x='cell number',y='bacteria RPM bins')+
        facet_grid(type~group,scale='free')+theme_classic()
    print(p_reads);print(p_RPM)
}

if (!args[8]%in%c('',NA)){
    species <- strsplit(args[8],',')[[1]]
    vars <- paste0(substr(args[7],1,1),'_',gsub('[a-z]* ','.',species),'Reads')
    sapply(species,summarise_reads)
    sapply(vars,scatter_plot)
    pdf(str_c(args[6],'/',sample,'.',k2db,'.hist.pdf'),width=9,height=9)
        plot_bin(vars)
    dev.off()
    p<-group.sum%>%filter(group!='no_tag')%>%
        pivot_longer(cols=c(vars),names_to='species',values_to='reads')%>%
        mutate(rpm=1e6*reads/TotalReads)%>%
        ggplot(aes(x=group,y=rpm,fill=species))+geom_col(position='dodge')+
        geom_text(aes(label=round(rpm,2)),position=position_dodge(0.9),vjust=-0)+labs(title=sample)
    ggsave(str_c(args[6],'/',sample,'.',k2db,'.RPM.pdf'),width=6,height=3)
}
k2out%>%filter(group!='no_tag')%>%      # ,!is.na(txID)
    select(c('corrected_barcode',ends_with('Reads'),ends_with('RPM')))%>%
    select(-kraken2Reads)%>%    # ,-totalReads,-mappedReads
    distinct()%>%write_tsv(str_c(args[6],'/',sample,'.',k2db,'.lg2a1.xls.gz'))
group.sum%>%write_tsv(str_c(args[6],'/',sample,'.',k2db,'.group.xls'))
scatter_plot('MicrobeReads')

p<-group.sum%>%filter(group!='no_tag')%>%ggplot(aes(x=group,y=Barcodes_counts))+
    geom_col()+geom_text(aes(label=Barcodes_counts),vjust=-0)+labs(title=sample)
ggsave(str_c(args[6],'/',sample,'.',k2db,'.cellfiltered.pdf'),width=2,height=3)
