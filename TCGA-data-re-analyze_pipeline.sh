#!/bin/bash
#SBATCH --mem-per-cpu=5000
#SBATCH --ntasks=4
source ~/software/miniforge3/bin/activate default
source scripts/bioFunction.sh.bak
# set -e
START=$(date +%s.%N)

# parameters:
    threads=4; host=hg38; peFq1suffix="_R1.fq.gz";seFqsuffix=".fq.gz"
    kk2DBlist=(~/databases/${host}_ABVF)
    # 后来陆续增加了一些需要看coverage的物种
    # microbeIndexes=(bacteriaGenomeIndex/{Alistipes_onderdonkii,Anaerotignum_propionicum,Bacteroides_fragilis,Cupriavidus_metallidurans,Curtobacterium_flaccumfaciens,Enterococcus_faecium,Fusobacterium_nucleatum,Methylobacterium_sp._C1,Mycobacterium_canettii,Staphylococcus_simulans,Streptococcus_anginosus})
    # microbeGffs=()
    # microbeIndexes=(bacteriaGenomeIndex/{Malassezia_restricta,Candida_albicans,Human_immunodeficiency_virus_1,Mus_musculus_mobilized_endogenous_polytropic_provirus,Tequatrovirus_gee4507,Woodchuck_hepatitis_virus,Aspergillus_oryzae,Botrytis_cinerea,Puccinia_triticina,Methylobacterium_aquaticum,Methylobacterium_brachiatum,Methylobacterium_currus,Methylobacterium_radiotolerans})
    # microbeGffs=()
    # microbeIndexes=(~/project/refSeq/bt2Index/{Fusobacterium_nucleatum,Faecalibacterium_cf._prausnitzii_KLE1255,Cutibacterium_acnes})
    # microbeGffs=(~/project/refSeq/fnagff/{GCF_001457555.1_NCTC10562_genomic.gff,GCF_000166035.1_ASM16603v1_genomic.gff,GCF_000376705.1_ASM37670v1_genomic.gff})
    microbeIndexes=(~/project/refSeq/bt2Index/{Staphylococcus_aureus,Klebsiella_pneumoniae,Prevotella_melaninogenica,Porphyromonas_asaccharolytica,Enterococcus_faecalis,Alphapapillomavirus_9,Streptococcus_mitis,Streptococcus_anginosus,Prevotella_intermedia,Streptococcus_pneumoniae,Human_gammaherpesvirus_4,Prevotella_jejuni,Veillonella_rogosae,Fusobacterium_nucleatum,Parvimonas_micra,Campylobacter_showae,Hungatella_hathewayi,Gemella_morbillorum,Campylobacter_curvus,Veillonella_nakazawae,Bacteroides_fragilis,Helicobacter_pylori,Alistipes_finegoldii,Veillonella_parvula,Veillonella_dispar,Phocaeicola_dorei})
    microbeGffs=(~/project/refSeq/fnagff/{GCF_000013425.1_ASM1342v1_genomic.gff,GCF_000240185.1_ASM24018v2_genomic.gff,GCF_000144405.1_ASM14440v1_genomic.gff,GCF_000212375.1_ASM21237v1_genomic.gff,GCF_000393015.1_Ente_faec_T5_V1_genomic.gff,GCF_000863945.3_ViralProj15505_genomic.gff,GCF_000960005.1_ASM96000v1_genomic.gff,GCF_001412635.1_ASM141263v1_genomic.gff,GCF_001953955.1_ASM195395v1_genomic.gff,GCF_002076835.1_ASM207683v1_genomic.gff,GCF_002402265.1_ASM240226v1_genomic.gff,GCF_002849795.1_ASM284979v1_genomic.gff,GCF_002959775.1_ASM295977v1_genomic.gff,GCF_003019295.1_ASM301929v1_genomic.gff,GCF_003454775.1_ASM345477v1_genomic.gff,GCF_004803815.1_ASM480381v1_genomic.gff,GCF_009721605.1_ASM972160v1_genomic.gff,GCF_009730315.1_ASM973031v1_genomic.gff,GCF_013372125.1_ASM1337212v1_genomic.gff,GCF_013393365.1_ASM1339336v1_genomic.gff,GCF_016889925.1_ASM1688992v1_genomic.gff,GCF_017821535.1_ASM1782153v1_genomic.gff,GCF_027677115.1_ASM2767711v1_genomic.gff,GCF_900186885.1_48903_D01_genomic.gff,GCF_900637515.1_51184_A01_genomic.gff,GCF_902387545.1_UHGG_MGYG-HGUT-02478_genomic.gff})

    picard=~/software/picard.jar
    hostIndex=~/project/refSeq/bt2Index/$host
    genome=~/project/refSeq/UCSC/$host.fa
    gtf=~/project/refSeq/UCSC/$host.ensGene.gtf
    if [ $host == 'mm39' ]; then gsize=mm; orgDb=org.Mm.eg.db;organism=mmu
    elif [ $host == 'hg38' ]; then gsize=hs; orgDb=org.Hs.eg.db;organism=hsa
    fi
# PIPELINE
    MAPHOST $hostIndex $threads bam $peFq1suffix $seFqsuffix
    [ -f sample.rename.group.csv ] || \
        awk 'BEGIN{FS="\t";OFS=","}{print $1,$1,$1}' merge/fq.map.xls | \
        sed '1csample,newname,group' >sample.rename.group.csv
    CLASSIFY "${kk2DBlist[*]}" $threads
    COVERAGE_INTERSECT "${microbeIndexes[*]}" "${microbeGffs[*]}" "${kk2DBlist[*]}" $threads
    GENERATE_RESULT "xls kraken2 coverage" $threads $genome $gtf $orgDb $organism
END=$(date +%s.%N)
date
echo "$END - $START" | bc
