args <- c('sample.rename.group.csv',
        'merge/fq.fastp.xls',
        'merge/kraken2.hg38_ABVF.taxon.xls',
        'merge/kraken2.hg38_ABVF.xls',
        'result/microbes/hg38_ABVF')
args <- commandArgs(trailingOnly = TRUE);args
sample.group.file <- args[[1]]
fq.totalReads.file <- args[[2]]
kk2taxon.file <- args[[3]]
kk2report.file <- args[[4]]
out.dir <- args[[5]]

library(tidyverse)
library(ggpubr)
library(ggpmisc)
library(stringr)
library(scales)
library(ggsci)
library(pheatmap)

# 读取样品分组信息
sample.group <- read_csv(sample.group.file)
# totalReads数值来自fq.fastp或者fq.map文件，自己测序文件必然有fastp，TCGA下载的bam文件只有map结果
# 注意：双端reads按对计数，与kraken2保持一致！
if (file.exists(fq.totalReads.file)){
    fq.totalReads <- read_tsv(fq.totalReads.file)%>%transmute(sample,totalReads=ReadsAfterFilt)
}else{
    fq.totalReads <- read_tsv(str_replace(fq.totalReads.file,'fastp','map'))%>%transmute(sample,totalReads=`total reads`)
    }
kk2taxon <- read_tsv(kk2taxon.file)%>%mutate(DBMinimizers=NULL,DBDistinctMinimizers=NULL); kk2taxon[is.na(kk2taxon)] <- "None"
kk2report <- read_tsv(kk2report.file) %>% filter(kraken2Minimizers>0 | kraken2ReadsCounts>0)    # unclassified 没有minimizer但有readsCounts
kk2dbName <- str_split(kk2taxon.file,"\\.")[[1]][2]

dir.create(out.dir,recursive = T)
setwd(out.dir)
a <- sample.group%>%left_join(fq.totalReads)%>%left_join(kk2report)%>%left_join(kk2taxon)%>%
    mutate(sample=factor(newname,levels=unique(newname)),newname=NULL,group=factor(group,levels=unique(group)))
# 计算Microbe Reads、Microbe minimizers、Microbe distinctminimizers
a <- a%>%group_by(sample)%>%filter(ScientificName%in%c("Fungi","Bacteria","Archaea","Viruses"))%>%
    summarise(
        kraken2MicrobeReads=sum(kraken2ReadsCounts),
        kraken2MicrobeMinimizers=sum(kraken2Minimizers),
        kraken2MicrobeDistinctMinimizers=sum(kraken2DistinctMinimizers)
    )%>%right_join(a)
tsv_col <- c('sample','group','totalReads','kraken2MicrobeReads','RankCode',
    'ScientificName','kraken2ReadsCounts','kraken2DistinctReadsCount','kraken2Minimizers','kraken2DistinctMinimizers',  # 
    'txID','Domain','Kindom','Phylum','Class','Order','Family','Genus','Species')
a[tsv_col] %>% write_tsv(basename(kk2report.file))
if (str_detect(kk2taxon.file,"_ABVF")){
    a[tsv_col]%>%filter(!RankCode%in%c('S','G'))%>%write_tsv(str_c('kraken2.',kk2dbName,'.DKPCOF.xls'))
    a[tsv_col]%>%filter(RankCode%in%c('R','U','S'))%>%write_tsv(str_c('kraken2.',kk2dbName,'.S.xls'))
    a[tsv_col]%>%filter(RankCode%in%c('R','U','G'))%>%write_tsv(str_c('kraken2.',kk2dbName,'.G.xls'))
    }

# 计算Bacteria，Fungi，Archaea，Viruses各自Reads比例
# 取RPM会导致RPM小于1的坐标经log变换后小于零
b <- a%>%filter(ScientificName%in%c("Fungi","Bacteria","Archaea","Viruses"))%>%
    transmute(group,sample,ScientificName,kraken2=kraken2ReadsCounts/totalReads)
    # pivot_longer(c(kraken2),names_to="Stat",values_to="ReadsPercent")

format1 <- function(plot){
    plot<-plot+
    scale_y_log10(
        limits = 10^c(-7,0),
        breaks = 10^c(-7:0),
        labels = trans_format("log10", math_format(10^.x)),
    )+
    # https://zhuanlan.zhihu.com/p/113551318
    # scale_y_continuous(
    #     trans = log10_trans(),
    #     breaks = trans_breaks("log10", function(x) 10^x),
    #     labels = trans_format("log10", math_format(10^.x)),
    #     limits = 10^c(-7,0),
    # )+
    labs(
        title="MicrobesReads / TotalReads",
        color="reads type",
        y="microbes reads / total reads"
        )+
    scale_color_d3(palette="category20")+theme_bw()+
    theme(
        panel.grid=element_line(color='gray97'),
        text=element_text(size=16,color='black'),
        axis.text=element_text(color='black'),
        plot.title=element_text(hjust = 0.5,size=18),
        axis.text.x=element_text(angle=90,hjust=1)
    )
}
m.group <- b%>%ggplot(aes(x=group,y=kraken2,color=ScientificName))+
    geom_boxplot(position=position_dodge2(width = 0.8))+
    geom_jitter(position=position_dodge2(width = 0.5),alpha=0.3)
format1(m.group)
width=max(4,0.5*length(levels(factor(b$group))))
ggsave("boxplot.MicrobesReads.group.pdf",height=4,width=width,limitsize=F)

m.sample <- b%>%ggplot(aes(x=sample,y=kraken2,color=ScientificName))+
    geom_boxplot(position="dodge")
format1(m.sample)
width=max(4,0.5*length(levels(factor(b$sample))))
ggsave("boxplot.MicrobesReads.sample.pdf",height=4,width=width,limitsize=F)

if (str_detect(kk2taxon.file,"_1Bc")){quit(save='no')}

# 计算不同Phylum、Genus的RPM比例柱形图
plot_bar <- function(rk,group,pst,type){
    p<-a%>%filter(RankCode==rk,
        !ScientificName%in%c("root","Chordata","Mus","Homo"))%>%
        mutate(kraken2=10^6*kraken2ReadsCounts/totalReads)%>%
        group_by(.data[[group]],ScientificName) %>% summarize(kraken2=mean(kraken2))%>%
        mutate(n=min_rank(desc(kraken2)),name=ifelse(n>3,"Others",ScientificName))%>%
        group_by(.data[[group]],name) %>% summarize(kraken2=sum(kraken2))%>%
        ggplot(aes(x=.data[[group]],y=kraken2,fill=name))+geom_col(position=pst)+
        labs(title=str_c(rk,".",type),fill="reads type",y=type)+
        # scale_y_continuous(labels = ifelse(pst=='fill',percent_format(),waiver()))+       # 默认waiver会出错
        scale_fill_d3(palette="category20")+theme_bw()+
        theme(panel.grid=element_line(color='gray97'),
            text=element_text(size=16,color='black'),
            axis.text=element_text(color='black'),
            plot.title=element_text(hjust = 0.5,size=18),
            axis.text.x=element_text(angle=90,hjust=1))
    width=max(4,0.5*length(levels(factor(a[[group]]))))
    ggsave(str_c("barplot.",rk,".",type,".",group,".pdf"),height=4,width=width,limitsize=F)
}
for(rk in c("P","G")){
    for(group in c('group','sample')){
        plot_bar(rk,group,'fill','abundance')
        plot_bar(rk,group,'stack','RPM')
    }
}

# Bacterial abundance Heatmap
if (str_detect(kk2taxon.file,"_21Bc")){RK=c("S1","G")}else{RK=c("S","G")}
for(rk in RK){
    t <- a %>% filter(RankCode==rk,
        !ScientificName%in%c("root","Homo sapiens","Mus musculus","Mus","Homo"),kraken2ReadsCounts>0)
    if (str_detect(kk2taxon.file,"_1Bc")){break}
    else if (str_detect(kk2taxon.file,"_21Bc")){
        c<-t%>%group_by(group,ScientificName)%>%summarize(kraken2ReadsCounts=mean(kraken2ReadsCounts),
                kraken2MicrobeReads=mean(kraken2MicrobeReads))%>%
            mutate(BcRatio=kraken2ReadsCounts/sum(kraken2ReadsCounts))%>%
            transmute(ScientificName,group,BcRatio)%>%
            pivot_wider(names_from=group,values_from=BcRatio,values_fill=0)
        c<-read_tsv(str_c("~/script/species_percentage.",rk,".xls"))%>%
            left_join(c,by=c("ScientificName"="ScientificName"))
        c<-c%>%filter(!is.na(ScientificName))%>%
            data.frame(row.names=.$ScientificName)%>%mutate(ScientificName=NULL)
        c[is.na(c)]<-0;c<-log10(c);c[c==-Inf]<--7

        d<-t%>%group_by(sample)%>%mutate(BcRatio=kraken2ReadsCounts/sum(kraken2ReadsCounts))%>%
            transmute(ScientificName,sample,BcRatio)%>%
            pivot_wider(names_from=sample,values_from=BcRatio,values_fill=0)
        d<-read_tsv(str_c("~/script/species_percentage.",rk,".xls"))%>%
            left_join(d,by=c("ScientificName"="ScientificName"))
        d<-d%>%filter(!is.na(ScientificName))%>%
            data.frame(row.names=.$ScientificName,check.names=F)%>%mutate(ScientificName=NULL)
        d[is.na(d)]<-0;d<-log10(d);d[d==-Inf]<--7
        # 按sample统计的group列注释，注意：由tibble转为dataframe会导致不规范列名发生改变
        anot_col_d<-a%>%filter(RankCode==rk,!ScientificName%in%c("root","Homo sapiens","Mus musculus"))%>%
            transmute(group,sample)%>%distinct()%>%data.frame(row.names=.$sample)%>%mutate(sample=NULL)
        anot_col_d<-rbind(data.frame(row.names='percentage',group='percentage'),anot_col_d)

        bk<-c(seq(-7,-2,by=0.01),seq(-1.99,0,by=0.01))
        pheatmap(c,cluster_cols=F,cluster_rows=F,angle_col="90",scale = "none",breaks=bk,
                color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/7*5),
                colorRampPalette(colors = c("white","orange","red"))(length(bk)/7*2)),
                main=str_c("21Bc/sum(21Bc), level=",rk),cellwidth=20,cellheight=12,
                filename=str_c('pheatmap.',rk,'.abundance.group.pdf'))
        pheatmap(d,annotation_col=anot_col_d,angle_col="90",scale = "none",breaks=bk,
                cluster_cols=F,cluster_rows=F,gaps_col=cumsum(rle(as.vector(anot_col_d$group))$lengths),
                color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/7*5),
                colorRampPalette(colors = c("white","orange","red"))(length(bk)/7*2)),
                main=str_c("21Bc/sum(21Bc), level=",rk),cellwidth=20,cellheight=12,
                filename=str_c('pheatmap.',rk,'.abundance.sample.pdf'))
    }
    else {      # if (str_detect(kk2taxon.file,"_ABVF"))
        # 按group统计
        c <- t%>%group_by(group,ScientificName)%>%summarize(kraken2ReadsCounts=mean(kraken2ReadsCounts),
                kraken2MicrobeReads=mean(kraken2MicrobeReads))%>%
            transmute(ScientificName,group,BcRatio=kraken2ReadsCounts/kraken2MicrobeReads)%>%
            pivot_wider(names_from=group,values_from=BcRatio,values_fill=0)
        c<-c%>%filter(!is.na(ScientificName))%>%
            data.frame(row.names=.$ScientificName,check.names=F)%>%mutate(ScientificName=NULL)
        c[is.na(c)]<-0;c<-log10(c);c[c==-Inf]<--7
        c<-c[apply(c,1,max)>-4,,drop=F]
        # 按sample统计
        d <- t%>%transmute(ScientificName,sample,BcRatio=kraken2ReadsCounts/kraken2MicrobeReads)%>%
            pivot_wider(names_from=sample,values_from=BcRatio,values_fill=0)
        d<-d%>%filter(!is.na(ScientificName))%>%
            data.frame(row.names=.$ScientificName,check.names=F)%>%mutate(ScientificName=NULL)
        d[is.na(d)]<-0;d<-log10(d);d[d==-Inf]<--7
        d<-d[apply(d,1,max)>-4,,drop=F]

        # 行注释，上面被筛选过，需要调整
        anot_row_cd<-a%>%filter(RankCode==rk,ScientificName%in%rownames(c),kraken2ReadsCounts>0)%>%
            transmute(Phylum,ScientificName)%>%distinct()%>%
            data.frame(row.names=.$ScientificName)%>%mutate(ScientificName=NULL)
        # 按sample统计的group列注释，注意：由tibble转为dataframe会导致不规范列名发生改变
        anot_col_d<-a%>%filter(RankCode==rk,!ScientificName%in%c("root","Homo sapiens","Mus musculus"),
            kraken2ReadsCounts>0)%>%
            transmute(group,sample)%>%distinct()%>%data.frame(row.names=.$sample)%>%mutate(sample=NULL)
        # rownames(anot_row_cd)==rownames(d); rownames(anot_col_d)==colnames(d)
        pheatmap(c,annotation_row=subset(anot_row_cd,rownames(anot_row_cd)%in%rownames(c)),
            cluster_cols=ifelse(length(ncol(c))>1,F,F),cluster_rows=ifelse(length(nrow(c))>1,T,F),angle_col="90",
            main=str_c("log10 Bacteria / Microbe, Max>-4, level=",rk),cellwidth=20,cellheight=12,
            filename=str_c('pheatmap.',rk,'.abundance.group.pdf'))
        pheatmap(d,annotation_row=subset(
            anot_row_cd,rownames(anot_row_cd)%in%rownames(d)),annotation_col=anot_col_d,
            cluster_cols=ifelse(length(ncol(d))>1,F,F),cluster_rows=ifelse(length(nrow(d))>1,T,F),angle_col="90",
            main=str_c("log10 Bacteria / Microbe, Max>-4, level=",rk),cellwidth=20,cellheight=12,
            gaps_col=cumsum(rle(as.vector(anot_col_d$group))$lengths),
            filename=str_c('pheatmap.',rk,'.abundance.sample.pdf'))

        c <- t%>%group_by(group,ScientificName)%>%summarize(kraken2ReadsCounts=mean(kraken2ReadsCounts),
                kraken2MicrobeReads=mean(kraken2MicrobeReads),totalReads=mean(totalReads))%>%
            transmute(ScientificName,group,BcRPM=10^6*kraken2ReadsCounts/totalReads)%>%
            pivot_wider(names_from=group,values_from=BcRPM,values_fill=0)
        c<-c%>%filter(!is.na(ScientificName))%>%
            data.frame(row.names=.$ScientificName,check.names=F)%>%mutate(ScientificName=NULL)
        c[is.na(c)]<-0;c<-log10(c);c[c==-Inf]<--1
        c<-c[apply(c,1,max)>0,,drop=F]
        # 按sample统计
        d <- t%>%transmute(ScientificName,sample,BcRPM=10^6*kraken2ReadsCounts/totalReads)%>%
            pivot_wider(names_from=sample,values_from=BcRPM,values_fill=0)
        d<-d%>%filter(!is.na(ScientificName))%>%
            data.frame(row.names=.$ScientificName,check.names=F)%>%mutate(ScientificName=NULL)
        d[is.na(d)]<-0;d<-log10(d);d[d==-Inf]<--1
        d<-d[apply(d,1,max)>0,,drop=F]
        pheatmap(c,annotation_row=subset(anot_row_cd,rownames(anot_row_cd)%in%rownames(c)),
            cluster_cols=ifelse(length(ncol(c))>1,F,F),cluster_rows=ifelse(length(nrow(c))>1,T,F),angle_col="90",
            main=str_c("log10 Bacteria RPM, Max>0, level=",rk),cellwidth=20,cellheight=12,
            filename=str_c('pheatmap.',rk,'.RPM.group.pdf'))
        pheatmap(d,annotation_row=subset(
            anot_row_cd,rownames(anot_row_cd)%in%rownames(d)),annotation_col=anot_col_d,
            cluster_cols=ifelse(length(ncol(d))>1,F,F),cluster_rows=ifelse(length(nrow(d))>1,T,F),angle_col="90",
            main=str_c("log10 Bacteria RPM, Max>0, level=",rk),cellwidth=20,cellheight=12,
            gaps_col=cumsum(rle(as.vector(anot_col_d$group))$lengths),
            filename=str_c('pheatmap.',rk,'.RPM.sample.pdf'))
    }
}

# 计算Genus层的c("Fungi","Bacteria","Archaea","Viruses")reads数（展敏去污以后）
reagentGenuslist <-c(
    "Aquabacterium",
    "Asticcacaulis",
    "Aurantimonas",
    "Beijerinckia",
    "Bosea",
    "Bradyrhizobium",
    "Brevundimonas",
    "Caulobacter",
    "Craurococcus",
    "Devosia",
    "Hoeflea",
    "Mesorhizobium",
    "Novosphingobium",
    "Paracoccus",
    "Pedomicrobium",
    "Phyllobacterium",
    "Rhizobium",
    "Sphingobium",
    "Sphingopyxis",
    "Acidovorax",
    "Azoarcus",
    "Azospira",
    "Curvibacter",
    "Duganella",
    "Herbaspirillum",
    "Janthinobacterium",
    "Limnobacter",
    "Massilia",
    "Methylophilus",
    "Methyloversatilis",
    "Pelomonas",
    "Polaromonas",
    "Ralstonia",
    "Schlegelella",
    "Sulfuritalea",
    "Undibacterium",
    "Variovorax",
    "Enhydrobacter",
    "Nevskia",
    "Pseudoxanthomonas",
    "Psychrobacter",
    "Xanthomonas",
    "Aeromicrobium",
    "Arthrobacter",
    "Beutenbergia",
    "Brevibacterium",
    "Curtobacterium",
    "Geodermatophilus",
    "Janibacter",
    "Microlunatus",
    "Patulibacter",
    "Rhodococcus",
    "Brevibacillus",
    "Brochothrix",
    "Paenibacillus",
    "Chryseobacterium",
    "Dyadobacter",
    "Flavobacterium",
    "Hydrotalea",
    "Niastella",
    "Olivibacter",
    "Pedobacter",
    "Deinococcus",
    "Afipia",
    "Comamonas",
    "Cupriavidus",
    "Delftia",
    "Leptothrix",
    "Oxalobacter",
    "Acinetobacter",
    "Stenotrophomonas",
    "Dietzia",
    "Kocuria",
    "Microbacterium",
    "Micrococcus",
    "Propionibacterium",
    "Tsukamurella",
    "Bacillus",
    "Mycobacterium",
    "Methylobacterium",
    "Ochrobactrum",
    "Roseomonas",
    "Sphingomonas",
    "Burkholderia",
    "Kingella",
    "Enterobacter",
    "Escherichia",
    "Pseudomonas",
    "Corynebacterium",
    "Abiotrophia",
    "Facklamia",
    "Streptococcus",
    "Wautersiella",
    "Staphylococcus"
)
if (str_detect(kk2taxon.file,"_ABVF")){
    dir.create("removed_contamination")
    setwd("removed_contamination")
    # 执行去污，并重新计算去污后的microbereads
    a_rc<-a%>%filter(RankCode=="G",kraken2ReadsCounts/totalReads*1000000>0.02,!Genus%in%reagentGenuslist)
    a_rc <- a_rc%>%group_by(sample)%>%filter(Kindom=="Fungi" | Domain%in%c("Bacteria","Archaea","Viruses"))%>%
        summarise(kraken2MicrobeReads_rc=sum(kraken2ReadsCounts))%>%right_join(a_rc)
    rc_col <- c('sample','group','totalReads','kraken2MicrobeReads','kraken2MicrobeReads_rc','RankCode',
        'ScientificName','kraken2ReadsCounts','txID','Domain','Kindom','Phylum','Class','Order','Family','Genus')
    a_rc[rc_col]%>%write_tsv(str_c('kraken2.',kk2dbName,'.removed_contamination.G.xls'))
    # 分别计算A,B,V,F的reads数量，准备boxplot
    abv<-a_rc%>%filter(Domain%in%c("Bacteria","Archaea","Viruses"))%>%
        group_by(group,sample,totalReads,Domain)%>%
        summarize(kraken2ReadsCounts=sum(kraken2ReadsCounts))%>%mutate(ScientificName=Domain)
    fungi<-a_rc%>%filter(Kindom=="Fungi")%>%group_by(group,sample,totalReads,Kindom)%>%
        summarize(kraken2ReadsCounts=sum(kraken2ReadsCounts))%>%mutate(ScientificName=Kindom)
    b<-full_join(abv,fungi)%>%filter(ScientificName%in%c("Fungi","Bacteria","Archaea","Viruses"))%>%
        transmute(group,sample,ScientificName,kraken2=kraken2ReadsCounts/totalReads)
    m.group <- b%>%ggplot(aes(x=group,y=kraken2,color=ScientificName))+geom_boxplot()
    format1(m.group)
    width=max(7,0.5*length(levels(factor(b$group))))
    ggsave("boxplot.MicrobesReads.group.pdf",width=width,limitsize=F)
    m.sample <- b%>%ggplot(aes(x=sample,y=kraken2,color=ScientificName))+geom_boxplot()
    format1(m.sample)
    width=max(7,0.5*length(levels(factor(b$sample))))
    ggsave("boxplot.MicrobesReads.sample.pdf",width=width,limitsize=F)

    rk="G"; a<-a_rc
    p.group<-a%>%filter(RankCode==rk,
        !ScientificName%in%c("root","Chordata","Mus","Homo"))%>%
        mutate(kraken2=kraken2ReadsCounts/totalReads)%>%
        group_by(group,ScientificName) %>% summarize(kraken2=mean(kraken2))%>%
        mutate(n=min_rank(desc(kraken2)),name=ifelse(n>3,"Others",ScientificName))%>%
        group_by(group,name) %>% summarize(kraken2=sum(kraken2))%>%
        ggplot(aes(x=group,y=kraken2,fill=name))+geom_col(position="fill")+
        labs(title=str_c(rk,".percentage"),fill="reads type",y="percentage")+
        scale_y_continuous(labels = percent_format())+
        scale_fill_d3(palette="category20")+theme_bw()+
        theme(panel.grid=element_line(color='gray97'),
            text=element_text(size=16),
            plot.title=element_text(hjust = 0.5,size=24),
            axis.text.x=element_text(angle=90,hjust=1))
    width=max(7,0.5*length(levels(factor(a$group))))
    ggsave(str_c("barplot.",rk,".percentage.group.pdf"),width=width,limitsize=F)

    p.sample<-a%>%filter(RankCode==rk,
        !ScientificName%in%c("root","Chordata","Mus","Homo"))%>%
        mutate(kraken2=kraken2ReadsCounts/totalReads)%>%
        group_by(sample,ScientificName) %>% summarize(kraken2=mean(kraken2))%>%
        mutate(n=min_rank(desc(kraken2)),name=ifelse(n>3,"Others",ScientificName))%>%
        group_by(sample,name) %>% summarize(kraken2=sum(kraken2))%>%
        ggplot(aes(x=sample,y=kraken2,fill=name))+geom_col(position="fill")+
        labs(title=str_c(rk,".percentage"),fill="reads type",y="percentage")+
        scale_y_continuous(labels = percent_format())+
        scale_fill_d3(palette="category20")+theme_bw()+
        theme(panel.grid=element_line(color='gray97'),
            text=element_text(size=16),
            plot.title=element_text(hjust = 0.5,size=24),
            axis.text.x=element_text(angle=90,hjust=1))
    width=max(7,0.5*length(levels(factor(a$sample))))
    ggsave(str_c("barplot.",rk,".percentage.sample.pdf"),width=width,limitsize=F)

    t <- a %>% filter(RankCode==rk,
        !ScientificName%in%c("root","Homo sapiens","Mus musculus","Mus","Homo"),kraken2ReadsCounts>0)
    # 按group统计
    c <- t%>%group_by(group,ScientificName)%>%summarize(kraken2ReadsCounts=mean(kraken2ReadsCounts),
            kraken2MicrobeReads=mean(kraken2MicrobeReads))%>%
        transmute(ScientificName,group,BcRatio=kraken2ReadsCounts/kraken2MicrobeReads)%>%
        pivot_wider(names_from=group,values_from=BcRatio,values_fill=0)
    c<-c%>%filter(!is.na(ScientificName))%>%
        data.frame(row.names=.$ScientificName,check.names=F)%>%mutate(ScientificName=NULL)
    c[is.na(c)]<-0;c<-log10(c);c[c==-Inf]<--7
    # 按sample统计
    d <- t%>%transmute(ScientificName,sample,BcRatio=kraken2ReadsCounts/kraken2MicrobeReads)%>%
        pivot_wider(names_from=sample,values_from=BcRatio,values_fill=0)
    d<-d%>%filter(!is.na(ScientificName))%>%
        data.frame(row.names=.$ScientificName,check.names=F)%>%mutate(ScientificName=NULL)
    d[is.na(d)]<-0;d<-log10(d);d[d==-Inf]<--7

    # 行注释，上面被筛选过，需要调整
    anot_row_cd<-a%>%filter(RankCode==rk,ScientificName%in%rownames(c),kraken2ReadsCounts>0)%>%
        transmute(Phylum,ScientificName)%>%distinct()%>%
        data.frame(row.names=.$ScientificName)%>%mutate(ScientificName=NULL)
    # 按sample统计的group列注释，注意：由tibble转为dataframe会导致不规范列名发生改变
    anot_col_d<-a%>%filter(RankCode==rk,!ScientificName%in%c("root","Homo sapiens","Mus musculus"),
        kraken2ReadsCounts>0)%>%
        transmute(group,sample)%>%distinct()%>%data.frame(row.names=.$sample)%>%mutate(sample=NULL)
    # rownames(anot_row_cd)==rownames(d); rownames(anot_col_d)==colnames(d)
    pheatmap(c,annotation_row=subset(anot_row_cd,rownames(anot_row_cd)%in%rownames(c)),
        cluster_cols=ifelse(length(ncol(c))>1,F,F),cluster_rows=T,angle_col="90",
        main=str_c("Bacteria / Microbe, Max>-4, level=",rk),cellwidth=20,cellheight=12,
        filename=str_c('pheatmap.',rk,'.abundance.group.pdf'))
    pheatmap(d,annotation_row=subset(
        anot_row_cd,rownames(anot_row_cd)%in%rownames(d)),annotation_col=anot_col_d,
        cluster_cols=ifelse(length(ncol(d))>1,F,F),cluster_rows=T,angle_col="90",
        main=str_c("Bacteria / Microbe, Max>-4, level=",rk),cellwidth=20,cellheight=12,
        gaps_col=cumsum(rle(as.vector(anot_col_d$group))$lengths),
        filename=str_c('pheatmap.',rk,'.abundance.sample.pdf'))
}
