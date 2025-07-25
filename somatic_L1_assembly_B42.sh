#!/bin/bash

export PATH=`pwd`/../mamba/bin:${PATH}

bcftools query -f "%CHROM\t%POS\t%FAM_N\t%ITYPE_N" B42tumor.vcf.gz  | grep "^chr" | egrep "solo|partnered"  | sort -k1,1V -k2,2n | uniq | awk '{print $1"\t"$2-5000"\t"($2+5000);}' > l1.pos
bedtools merge -i l1.pos > l1.pos.tmp
mv l1.pos.tmp l1.pos

while read CHR START END
do
    echo ${CHR} ${START}
    samtools view -F 3840 -b ../../hg38/B42tumor.bam ${CHR}:${START}-${END} > l1.bam
    samtools index l1.bam
    delly lr -g hg38.fa -d l1.dump.gz -o l1.bcf l1.bam
    bcftools query -i 'STRLEN(ALT)>STRLEN(REF)+500 && STRLEN(ALT)<STRLEN(REF)+9000' -f "%ID\n" l1.bcf > l1.svid
    if [ `cat l1.svid | wc -l` -eq 1 ]
    then
        SVID=`head -n 1 l1.svid`
        zcat l1.dump.gz | grep "^${SVID}" | cut -f 3 > l1.reads
        samtools view -h -N l1.reads B42tumor.bam ${CHR}:${START}-${END} | samtools fastq -TMM,ML - > l1.fq
        ./shasta-Linux-0.12.0 --input l1.fq --config Nanopore-May2022 --assemblyDirectory shasta.ass
        minimap2 -x map-ont -a -t 8 -y --secondary=no shasta.ass/Assembly.fasta l1.fq | samtools sort -m 3G -o assembly.bam -
        samtools index assembly.bam
        minimap2 -x map-ont -a -t 8 -y --secondary=no shasta.ass/Assembly.fasta LINE1.fa | samtools sort -m 3G -o line1.bam -
        samtools index line1.bam
        alfred bam2match -r shasta.ass/Assembly.fasta -o ${CHR}_${START}_${END}_line1.loc.in.assembly.gz line1.bam
        cp shasta.ass/Assembly.fasta ${CHR}_${START}_${END}_assembly.fasta
        cp assembly.bam ${CHR}_${START}_${END}_assembly.bam
        cp assembly.bam.bai ${CHR}_${START}_${END}_assembly.bam.bai
        cp line1.bam ${CHR}_${START}_${END}_line1.bam
        cp line1.bam.bai ${CHR}_${START}_${END}_line1.bam.bai
    fi
    rm -rf shasta.ass* 
done < l1.pos
