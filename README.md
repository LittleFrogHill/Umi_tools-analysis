# Umi_tools-analysis
analysis the reads with barcode

## 01.Umi_tools extract
using the raw data and put the Umi into the name of the reads

   for i in A006200178_153621_S1 \
    A006200178_153622_S2 \
    A006200178_153623_S3 \
    A006200178_153624_S4 \
    A006200178_153625_S5 \
    A006200178_153626_S6
    do
    umi_tools extract --bc-pattern=NNNNNNNNNNN --stdin=${i}_L002_R2_001.fastq.gz --read2-in=${i}_L002_R1_001.fastq.gz --stdout=${i}_add_barcode_R1.fastq.gz --read2-stdout
    umi_tools extract --bc-pattern=NNNNNNNNNNN --stdin=${i}_L002_R2_001.fastq.gz --read2-in=${i}_L002_R3_001.fastq.gz --stdout=${i}_add_barcode_R3.fastq.gz --read2-stdout
    done

## 02.clean data
using the specific adaptor seq 




bedtools intersect -a merge.snp.f.bi.vcf.gz -b /home/shangao/Data/EDTA/Ppr/Ppr/Ppr.FINAL.fa.mod.EDTA.TEanno.gff3 -v > merge.snp.f.bi.remove.repeat.vcf.gz
https://evomics.org/learning/population-and-speciation-genomics/2016-population-and-speciation-genomics/fileformats-vcftools-plink/#ex2.1
