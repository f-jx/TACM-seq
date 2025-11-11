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
    # 沙门:Salmonella_enterica GCF_000006945.2_ASM694v2_genomic.gff;肺克:Klebsiella_pneumoniae GCF_000240185.1_ASM24018v2_genomic.gff;
    # 金葡:Staphylococcus_aureus GCF_000013425.1_ASM1342v1_genomic.gff;Fuso:Fusobacterium_nucleatum GCF_001457555.1_NCTC10562_genomic.gff;
    # 枯草:Bacillus_subtilis GCF_000009045.1_ASM904v1_genomic.gff;大肠:Escherichia_coli GCF_000005845.2_ASM584v2_genomic.gff;脆拟:Bacteroides_fragilis GCF_016889925.1_ASM1688992v1_genomic.gff;
    microbeIndexes=(~/project/refSeq/bt2Index/{Moloney_murine_leukemia_virus,Bacteroides_fragilis,Faecalibacterium_prausnitzii,Roseburia_hominis,Veillonella_rogosae,Bifidobacterium_adolescentis,Fusobacterium_nucleatum,Limosilactobacillus_fermentum,Prevotella_copri,Escherichia_coli,Akkermansia_muciniphila,Candida_albicans,Clostridioides_difficile,Saccharomyces_cerevisiae,Methanobrevibacter_smithii,Salmonella_enterica,Enterococcus_faecalis,Clostridium_perfringens,21Akkermansia_muciniphila,21Bacteroides_fragilis,21Bifidobacterium_adolescentis,21Candida_albican,21Clostridioides_difficile,21Clostridium_perfringens,21Enterococcus_faecalis,21Escherichia_coli_B1109,21Escherichia_coli_b2207,21Escherichia_coli_B3008,21Escherichia_coli_B766,21Escherichia_coli_JM109,21Faecalibacterium_prausnitzii,21Fusobacterium_nucleatum,21Lactobacillus_fermentum,21Methanobrevibacter_smithii,21Prevotella_corporis,21Roseburia_hominis,21Saccharomyces_cerevisiae,21Salmonella_enterica,21Veillonella_rogosae})
    microbeGffs=(~/project/refSeq/fnagff/{GCF_000854185.1_ViralProj15030_genomic.gff,GCF_016889925.1_ASM1688992v1_genomic.gff,GCF_000154385.1_ASM15438v1_genomic.gff,GCF_902387955.1_UHGG_MGYG-HGUT-02517_genomic.gff,GCF_002959775.1_ASM295977v1_genomic.gff,GCF_003030905.1_ASM303090v1_genomic.gff,GCF_001457555.1_NCTC10562_genomic.gff,GCF_022819245.1_ASM2281924v1_genomic.gff,GCF_020735445.1_ASM2073544v1_genomic.gff,GCF_000005845.2_ASM584v2_genomic.gff,GCF_009731575.1_ASM973157v1_genomic.gff,GCF_000182965.3_ASM18296v3_genomic.gff,GCF_018885085.1_ASM1888508v1_genomic.gff,GCF_000146045.2_R64_genomic.gff,GCF_000016525.1_ASM1652v1_genomic.gff,GCF_000006945.2_ASM694v2_genomic.gff,GCF_000393015.1_Ente_faec_T5_V1_genomic.gff,GCF_020138775.1_ASM2013877v1_genomic.gff})

    picard=~/software/picard.jar
    hostIndex=~/project/refSeq/bt2Index/$host
    genome=~/project/refSeq/UCSC/$host.fa
    gtf=~/project/refSeq/UCSC/$host.knownGene.gtf
    if [ $host == 'mm39' ]; then gsize=mm; orgDb=org.Mm.eg.db;organism=mmu
    elif [ $host == 'hg38' ]; then gsize=hs; orgDb=org.Hs.eg.db;organism=hsa
    fi
# PIPELINE
    CELLRANGER atac/rna $host $threads
    MAPHOST $hostIndex $threads fq/bam $peFq1suffix $seFqsuffix
    RNA rsemIndex/$host bowtie2/star rsemMap/host $peFq1suffix $seFqsuffix $genome $gtf $threads
    GATK_CALL_SNP "bowtie2/host/*.bam" "snpOut/*" $genome $threads
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
    SUMMARY_SC_K2OUT atac/rna "kraken2/hg38_ABVF/*.k2out" $host $threads Species 'Salmonella enterica'
    GENERATE_RESULT "xls kraken2 coverage deseq2 atac diffbind enrich motif" $threads $genome $gtf $orgDb $organism
    rm fq/**/*.{gz,bam} atac_*/f2q30chrM/*.bam{,.bai} bowtie2/{microbes*[^e]/*,atac}/*.bam{,.bai}
END=$(date +%s.%N)
date
echo "$END - $START" | bc
