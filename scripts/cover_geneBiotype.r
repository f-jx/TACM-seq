args <- c('sample.rename.group.csv',
        'merge/fq.map.xls',
        'merge/merge.coverage.xls',
        'result/microbes/coverage')
args <- commandArgs(trailingOnly=TRUE);args
sample.group.file <- args[[1]]
fq.map.file <- args[[2]]
coverage.file <- args[[3]]
out.dir <- args[[4]]

# ???是否要把bases对测序深度做normalize，不用，因为计算coverage时，base和readscount的最小单位是1
# 一些高级用法：
# a%>%filter(!grepl("merge",sample),group==levels(factor(a$group))[[1]])
# com2<-com1%>%pivot_longer(names(com1)[-1:-4],names_to='stat',values_to='reads')
# if (length(levels(factor(a$group)))==0){next;quit(save='no')}
# # 删掉所有值都是na的列(kk2ReadsCounts有4个DB，有的DB不含有这个菌)
# if (length(which(apply(a,2,function(x)all(is.na(x)))))>0){a<-a[,-which(apply(a,2,function(x)all(is.na(x))))]}
# transmute(subproject,sample,species,
#     across(c(starts_with("kk2ReadsCounts"),bowtie2ReadsCount),function(x)(sprintf("%0.3f",10^6*x/RawTotalReads))))
#     names(k)<-gsub("bowtie2ReadsCount",str_replace(its,".*kk2Intersect","Intersect"),names(k))
library(tidyverse)
library(stringr)
library(scales)
library(ggrepel)
library(ggpmisc)
library(ggsci)

sample.group <- read_csv(sample.group.file)
fq.map <- read_tsv(fq.map.file)%>%transmute(sample,totalReads=`total reads`)
coverage <- read_tsv(coverage.file)
coverage <- sample.group%>%left_join(fq.map)%>%left_join(coverage)
names(coverage)
compare.reads <- coverage[c('sample','newname','group','totalReads','species','intersect','genomeReads')]%>%
    mutate(intersect=ifelse(intersect=='none','bowtie2',str_c('intersect.',intersect)))%>%
    pivot_wider(names_from=intersect,values_from=genomeReads,values_fill=NA)
for (kk2db in levels(factor(coverage$intersect))){
    if (kk2db=='none'){next}
    # kk2db='mm39_ABVF'
    species.txID <- read_tsv(str_c('merge/kraken2.',kk2db,'.taxon.xls'))%>%
        filter(RankCode==ifelse(str_detect(kk2db,'_21Bc'),'S1','S'))%>%transmute(txID,species=ScientificName)
    # 21菌的bowtie2结果和ABVF库17菌做交集，Ecoli合并成merge
    if(str_detect(kk2db,'_ABVF')){
        species.txID21 <- species.txID%>%filter(species%in%c('Bacteroides fragilis',
                                                            'Faecalibacterium prausnitzii',
                                                            'Roseburia hominis',
                                                            'Veillonella rogosae',
                                                            'Bifidobacterium adolescentis',
                                                            'Fusobacterium nucleatum',
                                                            'Limosilactobacillus fermentum',
                                                            'Prevotella copri',
                                                            'Escherichia coli',
                                                            'Akkermansia muciniphila',
                                                            'Candida albicans',
                                                            'Clostridioides difficile',
                                                            'Saccharomyces cerevisiae',
                                                            'Methanobrevibacter smithii',
                                                            'Salmonella enterica',
                                                            'Enterococcus faecalis',
                                                            'Clostridium perfringens'))
        species.txID21$species <- str_replace(species.txID21$species,'^Prevotella copri$','Prevotella corporis')
        species.txID21$species <- str_replace(species.txID21$species,'^Escherichia coli$','Ecoli merge')
        species.txID21$species <- str_replace(species.txID21$species,'^Candida albicans$','Candida albican')
        species.txID21$species <- str_replace(species.txID21$species,'^Limosilactobacillus fermentum$','Lactobacillus fermentum')
        species.txID21$species <- str_replace(species.txID21$species,'^','21')
        species.txID <- rbind(species.txID,species.txID21)
    }
    txID.kk2reads <- read_tsv(str_c('merge/kraken2.',kk2db,'.xls.gz'))%>%
        transmute(txID,sample,kraken2ReadsCounts)
    colnames(txID.kk2reads)[which(names(txID.kk2reads) == "kraken2ReadsCounts")] <- kk2db
    compare.reads <- compare.reads%>%left_join(species.txID)%>%
        left_join(txID.kk2reads)%>%mutate(txID=NULL)
}
dir.create(out.dir,recursive=TRUE)
setwd(out.dir)

compare.reads <- compare.reads%>%mutate(sample=newname,newname=NULL)
compare.reads %>% write_tsv('microbes.reads.compared.xls')

coverage <- coverage%>%mutate(sample=newname,newname=NULL)
coverage %>% write_tsv('microbes.coverage.xls')

# quit(save='no')
pdf('microbes.coverages.pdf',width=14)
for (sp in levels(factor(coverage$species))){
    # sp='21Salmonella enterica'; 
    print(sp)
    com1<-compare.reads%>%filter(species==sp,!grepl("16S",group))%>%
        pivot_longer(names(compare.reads)[-1:-4],names_to='type',values_to='reads',values_drop_na=T)%>%
        mutate(rpm=sprintf("%0.3f",10^6*reads/totalReads),
        reads=ifelse(reads==0,0.7,reads))%>%type_convert()%>%
        pivot_longer(c('reads','rpm'),names_to='stat',values_to='value')
p.reads<-com1%>%ggplot(aes(x=group,y=value,color=group))+
    geom_boxplot(position=position_dodge2(width = 0.8))+
    geom_jitter(position=position_dodge2(width = 0.5),alpha=0.3)+
    facet_grid(stat~type)+labs(title=str_c(sp," reads and rpm"))+
    # geom_text(aes(label=value),position=position_dodge(width=0.9),vjust=0)+
    # geom_text_repel(aes(label=ifelse(value>1,value,NA)),
    #     size=1.5, min.segment.length=unit(0.2, "lines"),
    #     position=position_dodge(width=0.9))+
    # scale_y_continuous(trans="sqrt",labels=scientific)+
    scale_y_continuous(trans = log10_trans(),
        breaks = trans_breaks("log10", function(x) 10^x),
        labels = trans_format("log10", math_format(10^.x)))+
    scale_color_d3(palette="category20")+
    theme_bw()+
    theme(panel.grid=element_line(color='gray97'),
        text=element_text(size=16),
        plot.title=element_text(hjust=0.5,size=20),
        axis.text.x=element_text(
            # size=min(16,700/nrow(com1)),
            angle=90,hjust=1))
print(p.reads)
# ggsave(str_c("test.pdf"),width=14)

    cov1<-coverage%>%filter(species==sp)%>%mutate(intersect=ifelse(intersect=='none','bowtie2',str_c('intersect.',intersect)))
p.coverage<-cov1%>%mutate(genomeCoverage=round(genomeCoverRatio_depth1,5))%>%
    ggplot(aes(x=genomeBases,y=genomeCoverage,color=group))+
    geom_point()+
    # geom_smooth(method="lm",se=FALSE)+
    scale_x_continuous(trans=log10_trans(),
        breaks=trans_breaks("log10", function(x) 10^x),
        labels=trans_format("log10", math_format(10^.x)))+
    scale_y_continuous(labels=percent_format(scale=100,
            accuracy=ifelse(min(na.rm=TRUE,cov1$genomeCoverage)>0.01,1,0.01)))+
    facet_wrap(~intersect)+
    labs(title=str_c(sp," genomeCoverage"))+
    scale_color_d3(palette="category20")+
    theme_bw()+
    theme(panel.grid=element_line(color='gray97'),
        text=element_text(size=16),
        plot.title=element_text(hjust=0.5,size=20))
print(p.coverage)
# ggsave(str_c("test.pdf"),width=14)

    biotype<-coverage%>%filter(species==sp)%>%
        mutate(intersect=ifelse(intersect=='none','bowtie2',str_c('intersect.',intersect)),
        r1.CDS=protein_codingBases,r2.rRNA=rRNABases,r3.others=othersBases)%>%
        pivot_longer(c(r1.CDS,r2.rRNA,r3.others),names_to="bioType",values_to="Bases")%>%
        group_by(group,intersect,bioType)%>%summarize(Bases=mean(Bases))
if(sum(!is.na(biotype$Bases))==0){next}
p.bt<-biotype%>%ggplot(aes(x=group,y=Bases,fill=bioType))+
    geom_col(position="fill")+
    scale_y_continuous(labels=percent_format(scale=100,accuracy=1))+
    # annotate(geom="table",x=4,y=3,label=bioTypeLengthPercent)+
    # geom_table(data,mapping=aes(x,y,label))
    facet_wrap(~intersect)+
    labs(title=str_c(sp," Bases bioType"))+
    # scale_fill_discrete(breaks=c("protein_coding","rRNA","others"),labels=c(pLP,rLP,oLP))+
    scale_fill_d3(palette="category20",)+
        # breaks=c("r1.CDS","r2.rRNA","r3.others"),labels=c(pLP,rLP,oLP))+
    theme_bw()+
    theme(panel.grid=element_line(color='gray97'),
        text=element_text(size=16),
        plot.title=element_text(hjust=0.5,size=20),
        axis.text.x=element_text(angle=90,hjust=1))
print(p.bt)
# ggsave(str_c("test.pdf"),width=14)

# # gene_biotype Coverage
# cov2<-coverage%>%mutate(
#     r4.genome=CoverPercent_depth1,
#     r1.CDS=(pcCoverLength/protein_codingLength),
#     r2.rRNA=(rRNACoverLength/rRNALength),
#     r3.others=(otsCoverLength/othersLength))%>%
#     group_by(group,intersect,species)%>%
#     summarize(across(c(r4.genome,r1.CDS,r2.rRNA,r3.others),mean))
# cov2%>%
#     pivot_longer(c(r4.genome,r1.CDS,r2.rRNA,r3.others),names_to="bioType",values_to="Coverage")

# cov2%>%mutate(
#         r1.CDS=r1.CDS/r4.genome,
#         r2.rRNA=r2.rRNA/r4.genome,
#         r3.others=r3.others/r4.genome)%>%
#     pivot_longer(c(r1.CDS,r2.rRNA,r3.others),names_to="bioType",values_to="RelativeCoverageToGenome")
}
dev.off()