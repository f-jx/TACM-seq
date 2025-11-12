* libraries.csv
```
fastqs,sample,library_type
/md01/nieyg/project/OSN_mm/data/negative/scRNA,negative,Gene Expression
/md01/nieyg/project/OSN_mm/data/negative/scATAC,negative,Chromatin Accessibility
```

* 01_Positive_cellranger.pbs
```
#PBS -N Positive
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=24
#PBS -l mem=64000m

source /public/home/nieyg/.bash_profile
cd /md01/nieyg/project/OSN_mm/data/cellranger
cellranger-arc count --id=Positive \
                       --reference=/md01/nieyg/ref/10X/refdata-cellranger-arc-mm10-2020-A-2.0.0 \
                       --libraries=/md01/nieyg/project/OSN_mm/data/positive/libraries.csv \
                       --localcores=24 \
                       --localmem=64
```

* 02_Negative_cellranger.pbs
```
#PBS -N Negative
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=24
#PBS -l mem=64000m

source /public/home/nieyg/.bash_profile
cd /md01/nieyg/project/OSN_mm/data/cellranger
cellranger-arc count --id=Negative \
                       --reference=/md01/nieyg/ref/10X/refdata-cellranger-arc-mm10-2020-A-2.0.0 \
                       --libraries=/md01/nieyg/project/OSN_mm/data/negative/libraries.csv \
                       --localcores=24 \
```
# filter RNA and ATAC bam files 

nohup sh filter_bam.sh &

# get fragments file from filtered bam file 
# 1. ATAC 
conda activate python38

newbam=atac_possorted_bam.dedup.bam
samtools index -@ 12 $newbam
sinto fragments -b $newbam -p 12 -f atac_fragments_filtered.tsv 
#sinto fragments -b filtered.bam -f fragments.tsv 

cat atac_fragments_filtered.tsv |sort -k1,1 -k2,2n | bgzip > atac_fragments_filtered.tsv.gz
tabix -b 2 -e 3 -p bed atac_fragments_filtered.tsv.gz

# 2. RNA 

# only keep the unique reads in bam file 

cp /md01/nieyg/project/honeybee/data/cellranger/Forager/Forager-NCBI-manually/outs/parse_cellsorted_bam.py .
cp /md01/nieyg/project/honeybee/data/cellranger/Forager/Forager-NCBI-manually/outs/run_dedup.sh .

nohup bash run_dedup.sh gex_possorted_bam &








