#!/bin/bash
#SBATCH --mem-per-cpu=5000
#SBATCH --ntasks=4
source ~/software/miniforge3/bin/activate default
source scripts/bioFunction.sh.bak
# set -e
START=$(date +%s.%N)

# parameters:
    threads=4; host=mm39; peFq1suffix="_R1.fq.gz";seFqsuffix=".fq.gz"
    kk2DBlist=(~/databases/${host}_{21Bc,ABVF})
    microbeIndexes=(~/project/refSeq/bt2Index/{21Akkermansia_muciniphila,21Bacteroides_fragilis,21Bifidobacterium_adolescentis,21Candida_albican,21Clostridioides_difficile,21Clostridium_perfringens,21Enterococcus_faecalis,21Escherichia_coli_B1109,21Escherichia_coli_b2207,21Escherichia_coli_B3008,21Escherichia_coli_B766,21Escherichia_coli_JM109,21Faecalibacterium_prausnitzii,21Fusobacterium_nucleatum,21Lactobacillus_fermentum,21Methanobrevibacter_smithii,21Prevotella_corporis,21Roseburia_hominis,21Saccharomyces_cerevisiae,21Salmonella_enterica,21Veillonella_rogosae})
    microbeGffs=()

    picard=~/software/picard.jar
    hostIndex=~/project/refSeq/bt2Index/$host
    genome=~/project/refSeq/UCSC/$host.fa
    gtf=~/project/refSeq/UCSC/$host.ensGene.gtf
    if [ $host == 'mm39' ]; then gsize=mm; orgDb=org.Mm.eg.db;organism=mmu
    elif [ $host == 'hg38' ]; then gsize=hs; orgDb=org.Hs.eg.db;organism=hsa
    fi
# PIPELINE
    MAPHOST $hostIndex $threads fq $peFq1suffix $seFqsuffix
    [ -f sample.rename.group.csv ] || \
        awk 'BEGIN{FS="\t";OFS=","}{print $1,$1,$1}' merge/fq.map.xls | \
        sed '1csample,newname,group' >sample.rename.group.csv
    CLASSIFY "${kk2DBlist[*]}" $threads
    COVERAGE_INTERSECT "${microbeIndexes[*]}" "${microbeGffs[*]}" "${kk2DBlist[*]}" $threads

    mkdir bowtie2/atac -p && ln -f bowtie2/host/* bowtie2/atac/
    ATAC filter $picard $gtf $gsize $threads
    for i in `ls atac_* -d|sed 's#atac_##'`;do
        [ -f diffbind_$i.csv ] || awk 'BEGIN{FS=OFS=","}NR==FNR{a[$1]=$2;b[$1]=$3}
            NR!=FNR && FNR==1{print $0}NR!=FNR && FNR>1{
                $2=b[$1];$1=a[$1];$3=substr($1,length($1));print $0}' \
            sample.rename.group.csv merge/diffbind_$i.csv > diffbind_$i.csv
    done
    GENERATE_RESULT "xls kraken2 coverage atac" $threads $genome $gtf $orgDb $organism
END=$(date +%s.%N)
date
echo "$END - $START" | bc
