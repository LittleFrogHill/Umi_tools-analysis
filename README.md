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

## 02.clean data and mapping
using the specific adaptor seq 
	   for i in A006200178_153621_S1 \
		A006200178_153622_S2 \
		A006200178_153623_S3 \
		A006200178_153624_S4 \
		A006200178_153625_S5 \
		A006200178_153626_S6
	do
		#A006200178_153621_S1 A006200178_153622_S2 \
		#~/Software/QC/TrimGalore-0.6.5/trim_galore -j 30 -q 30 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -a2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --fastqc --paired --output_dir clean_data ${i}_add_barcode_R1.fastq.gz ${i}_add_barcode_R3.fastq.gz 
		#bwa mem -t 40 -R "@RG\tID:${i}\tSM:${i}\tLB:${i}\tPL:Illumina" /home/shangao/Data/juicer/Ppr/test/review/blobtools/genome/Ppr.FINAL.fa clean_data/${i}_add_barcode_R1_val_1.fq.gz clean_data/${i}_add_barcode_R3_val_2.fq.gz > clean_data/${i}.sam
		#python /home/shangao/script/python/umi_tools_change_reads_title.py -s clean_data/${i}.sam -o clean_data/${i}.replace.sam
		#LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/:$LD_LIBRARY_PATH
		#export LD_LIBRARY_PATH
		#bowtie --threads 40 -v 2 -m 10 -a /home/shangao/Data/hifiasm_tell-seq/Ppr/Ppr -1 clean_data/${i}_add_barcode_R1_val_1.fq.gz -2 clean_data/${i}_add_barcode_R3_val_2.fq.gz  --sam > clean_data/${i}.sam
		#samtools view -bS clean_data/${i}.sam > clean_data/${i}.bam
		#samtools sort clean_data/${i}.bam -o clean_data/${i}.sort.bam
		#samtools index clean_data/${i}.sort.bam
		#umi_tools dedup -I clean_data/${i}.sort.bam --output-stats=deduplicated --paired -S clean_data/${i}.sort.de.bam --temp-dir=tmp
		#samtools flagstat clean_data/${i}.sort.de.bam > clean_data/${i}.sort.de.flag
	        #samtools flagstat clean_data/${i}.sort.bam > clean_data/${i}.sort.flag

		#samtools flagstat clean_data/${i}.sort.de.bam > ~/linda_mites/${i}.sort.de.flag
	        #samtools flagstat clean_data/${i}.sort.bam > ~/linda_mites/${i}.sort.flag
		#samtools index clean_data/${i}.sort.de.bam
		
	done

## 03. SNP calling
	   #echo "sh gatk.sh ../${i}.sort.de.bam ${i}.g.vcf ${i}.vcf" > clean_data/snp/script/$i.sh
		echo "sh filter.sh ${i}.vcf.gz ${i}.snp.vcf.gz ${i}.indel.vcf.gz ${i}.f.snp.vcf.gz ${i}.f.indel.vcf.gz ${i}.f.vcf.gz" > clean_data/snp/script/$i.f.sh
		#pigz -d clean_data/snp/$i.vcf.gz
		#bgzip -f clean_data/snp/$i.vcf
		#tabix -p vcf clean_data/snp/$i.vcf.gz

## 04. gatk.sh
	   cat gatk.sh 
### software
	gatk=/NVME/Software/popgen/gatk-4.1.9.0/gatk
### data
	ref=/home/shangao/Data/juicer/Ppr/test/review/blobtools/genome/Ppr.FINAL.fa
	bam=$1
	out1=$2
	out2=$3
### build dict for genome
	#/NVME/Software/popgen/gatk-4.1.9.0/gatk CreateSequenceDictionary -R Ppr.FINAL.fa -O Ppr.FINAL.dict
### call gvcf
	$gatk HaplotypeCaller \
		-R $ref \
		--emit-ref-confidence GVCF \
		-I $bam \
		-O $out1
	
### detect SNPs
	$gatk GenotypeGVCFs \
		-R $ref \
		-V $out1 \
		-O $out2

### compress
	bgzip -f $out2
	tabix -p vcf $out2.gz

## 05. cat filter.sh
	###software
	gatk=/NVME/Software/popgen/gatk-4.1.9.0/gatk
	###data
	ref=/home/shangao/Data/juicer/Ppr/test/review/blobtools/genome/Ppr.FINAL.fa
	vcf=$1
	snpvcf=$2
	indelvcf=$3
	filterSNP=$4
	filterINDEL=$5
	finalvcf=$6
### SelectVariants SNP
	$gatk SelectVariants \
		-select-type SNP \
		-V $vcf \
		-O $snpvcf
### SelectVariants INDEL
	$gatk SelectVariants \
	        -select-type INDEL \
	        -V $vcf \
	        -O $indelvcf
### filter SNP
	$gatk VariantFiltration \
		-V $snpvcf \
		--filter-expression "QD <2.0 || MQ <30.0 || FS >60.0 || SOR >5.0 || ReadPosRankSum < -8.0" \
		--filter-name "PASS" \
		-O $filterSNP
### filter INDEL
	$gatk VariantFiltration \
	        -V $indelvcf \
	        --filter-expression "QD <2.0 || FS >100.0 || SOR >5.0 || ReadPosRankSum < -8.0" \
	        --filter-name "PASS" \
	        -O $filterINDEL
### merge SNP INDEL
	$gatk MergeVcfs \
		-I $filterSNP \
		-I $filterINDEL \
		-O $finalvcf
### delete temp
	rm -f $snpvcf $indelvcf $filterSNP $filterINDEL

## 06.cat merged_gvcf.sh
	/NVME/Software/popgen/gatk-4.1.9.0/gatk CombineGVCFs -R /home/shangao/Data/juicer/Ppr/test/review/blobtools/genome/Ppr.FINAL.fa -V A006200178_153621_S1.g.vcf -V A006200178_153622_S2.g.vcf -V A006200178_153623_S3.g.vcf -V A006200178_153624_S4.g.vcf -V A006200178_153625_S5.g.vcf -V A006200178_153626_S6.g.vcf -O merged_gvcf/merge.g.vcf.gz
	/NVME/Software/popgen/gatk-4.1.9.0/gatk GenotypeGVCFs \
		-R /home/shangao/Data/juicer/Ppr/test/review/blobtools/genome/Ppr.FINAL.fa \
		-V merged_gvcf/merge.g.vcf.gz \
		-O merged_gvcf/merge.vcf.gz
	
	/NVME/Software/popgen/gatk-4.1.9.0/gatk SelectVariants \
		-select-type SNP \
		-V merged_gvcf/merge.vcf.gz \
		-O merged_gvcf/merge.snp.vcf.gz
	
	/NVME/Software/popgen/gatk-4.1.9.0/gatk SelectVariants \
		-select-type INDEL \
	        -V merged_gvcf/merge.vcf.gz \
	        -O merged_gvcf/merge.indel.vcf.gz
	
	/NVME/Software/popgen/gatk-4.1.9.0/gatk VariantFiltration \
		-V merged_gvcf/merge.snp.vcf.gz \
		--filter-expression "QD <2.0 || MQ <30.0 || FS >60.0 || SOR >5.0 || ReadPosRankSum < -8.0" \
		--filter-name "PASS" \
		-O merged_gvcf/merge.snp.f.vcf.gz
### filter INDEL
	/NVME/Software/popgen/gatk-4.1.9.0/gatk VariantFiltration \
	        -V merged_gvcf/merge.indel.vcf.gz \
	        --filter-expression "QD <2.0 || FS >100.0 || SOR >5.0 || ReadPosRankSum < -8.0" \
	        --filter-name "PASS" \
	        -O merged_gvcf/merge.indel.f.vcf.gz
### merge SNP INDEL
	/NVME/Software/popgen/gatk-4.1.9.0/gatk MergeVcfs \
		-I merged_gvcf/merge.snp.f.vcf.gz \
		-I merged_gvcf/merge.indel.f.vcf.gz \
		-O merged_gvcf/merge.f.vcf.gz
	
## 07.cat merged_gvcf/paired/extract_SNP.sh
	vcftools --gzvcf ../merge.snp.f.bi.vcf.gz \
		--recode-INFO-all \
		--maxDP 30 \
		--minDP 10 \
		--minQ 30 \
		--recode \
		--stdout \
		--maf 0.05 \
		--min-meanDP 20 \
		--max-missing 0.95 \
		--indv A006200178_153625_S5 \
		--indv A006200178_153626_S6 \
		--out s56.vcf > s56.vcf
	
	grep -v '|' s56.vcf > s56.removedshuxian.vcf
	
	bedtools intersect -a s56.removedshuxian.vcf -b /home/shangao/Data/EDTA/Ppr/Ppr/Ppr.FINAL.fa.mod.EDTA.TEanno.gff3 -v > s56.removedshuxian.afterbed.vcf
	
	python /home/shangao/script/python/vcf_filter-same.py -s s56.removedshuxian.afterbed.vcf -o s56.removedshuxian.afterbed.removesame.vcf

	python /home/shangao/script/python/vcf_filter-same.1.py -s s56.removedshuxian.afterbed.removesame.vcf -o s56.removedshuxian.afterbed.removesame.filterread.vcf

## 08.cat merged_gvcf/paired/pca/test.sh 
	vcftools --gzvcf ../../merge.snp.f.bi.vcf.gz --plink --out merge.snp.f.bi.remove.repeat.vcf
	/home/shangao/Software/plink --noweb --file merge.snp.f.bi.remove.repeat.vcf --make-bed --out merge.snp.f.bi.remove.repeat.vcf_bfile
	/home/shangao/Software/plink --threads 16 --bfile merge.snp.f.bi.remove.repeat.vcf_bfile --pca 3 --out merge.snp.f.bi.remove.repeat.vcf_pca3_bfile
	library(ggplot2)
	library(ggrepel)
	data<-read.table("merge.snp.f.bi.remove.repeat.vcf_pca3_bfile.eigenvec",header=T)
	> pdf('pca12.pdf')
	> ggplot(data,aes(x=pc1,y=pc2, colour = col,shape=col))+ geom_point()+ geom_text_repel(aes(x=pc1,y=pc2,label=name))+scale_x_continuous(limits=c(-1, 1))
	> dev.off()
	pdf
	  2
	> pdf('pca13.pdf')
	> ggplot(data,aes(x=pc1,y=pc3, colour = col,shape=col))+ geom_point()+ geom_text_repel(aes(x=pc1,y=pc3,label=name))+scale_x_continuous(limits=c(-1, 1))
	> dev.off()
	pdf
	  2
	> pdf('pca23.pdf')
	> ggplot(data,aes(x=pc2,y=pc3, colour = col,shape=col))+ geom_point()+ geom_text_repel(aes(x=pc2,y=pc3,label=name))+scale_x_continuous(limits=c(-1, 1))
	> dev.off()

## 09.cat merged_gvcf/paired/plot_dot/test.sh 
	grep 'chr1' ../s12.removedshuxian.afterbed.removesame.filterread.vcf|awk '{print$1,$2}' > s12.chr1.plot
	sed -i 's/chr1/1.5/g' s12.chr1.plot
	
	grep 'chr1' ../s34.removedshuxian.afterbed.removesame.filterread.vcf|awk '{print$1,$2}' > s34.chr1.plot
	sed -i 's/chr1/1/g' s34.chr1.plot
	
	grep 'chr1' ../s56.removedshuxian.afterbed.removesame.filterread.vcf|awk '{print$1,$2}' > s56.chr1.plot
	sed -i 's/chr1/0.5/g' s56.chr1.plot
	
	cat title s12.chr1.plot s34.chr1.plot s56.chr1.plot > chr1.plot
	
	library(ggplot2)
	data<-read.table("s12.plot",header=T)
	pdf('chr1.pdf',width=20,height=3)
	ggplot(data,aes(x=pos,y=col, colour = col))+ geom_point()+ scale_y_continuous(limits=c(0, 2))+ scale_x_continuous(limits=c(1, 31773567))+theme(panel.background = element_blank(),panel.border = element_blank(),legend.position="none",axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank(),plot.title = element_text(hjust = 0.5))+ggtitle("Chr1")+ylab("Position")
	dev.off()
	cat merged_gvcf/paired/plot_dot/test1.sh 
	python /home/shangao/script/python/calulate_SNp_num.py -s ../s56.removedshuxian.afterbed.removesame.filterread.vcf -o s56.slide

	pdf('All.pdf',width=20)
	ggplot(data,aes(x=pos,y=num, colour = chr,shape=as.factor(shape)))+ geom_point(size=3)+ scale_y_continuous(limits=c(0, 20),expand=c(0,0))+scale_x_continuous(expand=c(0,0))+theme(panel.background = element_blank(),panel.border = element_blank(),plot.title = element_text(hjust = 0.5),axis.line.x = element_line(color = "black"),axis.line.y = element_line(color = "black"))+ggtitle("Number of mutation")
	dev.off()
	
	
	
	bedtools intersect -a merge.snp.f.bi.vcf.gz -b /home/shangao/Data/EDTA/Ppr/Ppr/Ppr.FINAL.fa.mod.EDTA.TEanno.gff3 -v > merge.snp.f.bi.remove.repeat.vcf.gz
	https://evomics.org/learning/population-and-speciation-genomics/2016-population-and-speciation-genomics/fileformats-vcftools-plink/#ex2.1
