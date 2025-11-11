#!/bin/bash
#SBATCH --mem-per-cpu=5000
#SBATCH --ntasks=4
source ~/software/miniforge3/bin/activate default
source scripts/bioFunction.sh.bak
# set -e
START=$(date +%s.%N)

# parameters:
    threads=4; host=mm39; peFq1suffix="_R1.fq.gz";seFqsuffix=".fq.gz"
    kk2DBlist=(~/databases/${host}_ABVF)
    # microbeIndexes=(~/project/refSeq/bt2Index/Mycoplasma_pulmonis)
    # microbeGffs=()

    picard=~/software/picard.jar
    hostIndex=~/project/refSeq/bt2Index/$host
    genome=~/project/refSeq/UCSC/$host.fa
    gtf=~/project/refSeq/UCSC/$host.knownGene.gtf
    if [ $host == 'mm39' ]; then gsize=mm; orgDb=org.Mm.eg.db;organism=mmu
    elif [ $host == 'hg38' ]; then gsize=hs; orgDb=org.Hs.eg.db;organism=hsa
    fi

# BioFunction
    # CELLRANGER atac/rna hg38/mm39 $threads
    CELLRANGER(){
        local bam
        echo cellRanger flagstat bam2fq fastp
        CELLRANGER_COUNT $1 fq/sc bowtie2/host $2 $3
        [ $1 == 'rna' ] && bam="bowtie2/host/*/outs/possorted_genome_bam.bam"
        [ $1 == 'atac' ] && bam="bowtie2/host/*/outs/possorted_bam.bam"
        BAM_FLAGSTAT "$bam" merge/fq.map.xls $3
        # 需要对cellranger的bam进行二次filter
        SCBAMQC "bam" $3; BAM2FQ "${bam%%.bam}.filtered_out.bam" all "fq/unmap/raw/*" $3
        # 有时候cellranger的bam文件会出现reads长度为0的情况，需要对空行填充-，否则后续fastp会报错
        for fq in `ls fq/unmap/raw/*.gz`;do zcat $fq | sed 's#^$#-#' | gzip > $fq.tmp && mv $fq.tmp $fq; done
        FASTP fq/unmap/raw fq/unmap .1.gz .gz merge/unmapfq.fastp.xls $3
    }

# PIPELINE
    CELLRANGER atac $host $threads
    [ -f sample.rename.group.csv ] || \
        awk 'BEGIN{FS="\t";OFS=","}{print $1,$1,$1}' merge/fq.map.xls | \
        sed '1csample,newname,group' >sample.rename.group.csv
    CLASSIFY "${kk2DBlist[*]}" $threads
    # COVERAGE_INTERSECT "${microbeIndexes[*]}" "${microbeGffs[*]}" "${kk2DBlist[*]}" $threads

    SUMMARY_SC_K2OUT atac "kraken2/hg38_ABVF/*.k2out" $host $threads Species 'Mycoplasma pulmonis'
    GENERATE_RESULT "xls kraken2 coverage" $threads $genome $gtf $orgDb $organism
END=$(date +%s.%N)
date
echo "$END - $START" | bc
