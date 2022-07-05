# Dowloading sample dataset
echo -e "\n Downloading data... \n"
	
mkdir -p raw_data 
cd raw_data
	
wget 
https://zenodo.org/record/2582555/files/SLGFSK-N_231335_r1_chr5_12_17.fastq.gz
wget 
https://zenodo.org/record/2582555/files/SLGFSK-N_231335_r2_chr5_12_17.fastq.gz
wget 
https://zenodo.org/record/2582555/files/SLGFSK-T_231336_r1_chr5_12_17.fastq.gz
wget 
https://zenodo.org/record/2582555/files/SLGFSK-T_231336_r2_chr5_12_17.fastq.gz

#Downloading reference
echo -e "\n Downloading reference sequence... \n"
	
wget https://zenodo.org/record/2582555/files/hg19.chr5_12_17.fa.gz

#unzip reference 
gunzip hg19.chr5_12_17.fa.gz	

# Quality check
echo -e "\n Data Preprocessing... \n"

mkdir -p Fastqc_Reports  #create directory for the fastqc output

# creating list.txt ( a file containing the names of all the 4 files that 
have been dowloaded)
nano list.txt
  #type out the names of the samples and save buffer

#Qc on reads
for sample in `cat list.txt`
do
	fastqc ${sample}*.fastq.gz -o Fastqc_Reports
done

multiqc Fastqc_Reports -o Fastqc_Reports

cd ..

#dowloading the adapter sequence file
wget 
https://github.com/timflutre/trimmomatic/raw/master/adapters/TruSeq3-PE.fa 
-O TruSeq3-PE.fa

mkdir -p trimmed_reads

#Removing low quality sequences using Trimmomatic
for sample in `cat raw_data/list.txt`
do
       trimmomatic PE -threads 8 raw_data/${sample}_r1_chr5_12_17.fastq.gz 
raw_data/${sample}_r2_chr5_12_17.fastq.gz \
               trimmed_reads/${sample}_r1_paired.fq.gz 
trimmed_reads/${sample}_r1_unpaired.fq.gz \
               trimmed_reads/${sample}_r2_paired.fq.gz 
trimmed_reads/${sample}_r2_unpaired.fq.gz \
               ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:8:keepBothReads \
               LEADING:3 TRAILING:10 MINLEN:25
       
       fastqc  trimmed_reads/${sample}_r1_paired.fq.gz  
trimmed_reads/${sample}_r2_paired.fq.gz \
                 -o trimmed_reads/Fastqc_results
done 

multiqc  trimmed_reads/Fastqc_results  -o trimmed_reads/Fastqc_results

#Index reference file	
bwa index raw_data/hg19.chr5_12_17.fa

mkdir Mapping
   
#Perform alignment
bwa mem -R '@RG\tID:231335\tSM:Normal' raw_data/hg19.chr5_12_17.fa 
trimmed_reads/SLGFSK-N_231335_r1_paired.fq.gz \
      trimmed_reads/SLGFSK-N_231335_r2_paired.fq.gz > 
Mapping/SLGFSK-N_231335.sam

bwa mem -R '@RG\tID:231336\tSM:Tumor' raw_data/hg19.chr5_12_17.fa 
trimmed_reads/SLGFSK-T_231336_r1_paired.fq.gz \
       trimmed_reads/SLGFSK-T_231336_r2_paired.fq.gz > 
Mapping/SLGFSK-T_231336.sam

#Converting SAM to BAM, sorting and the output BAM.
for sample in `cat raw_data/list.txt`
do
        #Convert SAM to BAM and sort it 
        samtools view -@ 20 -S -b Mapping/${sample}.sam | samtools sort -@ 
32 > Mapping/${sample}.sorted.bam
        
        #Index BAM file
        samtools index Mapping/${sample}.sorted.bam
done

#Filtering Mapped reads
for sample in `cat raw_data/list.txt`
do
	#Filter BAM files
        samtools view -q 1 -f 0x2 -F 0x8 -b Mapping/${sample}.sorted.bam > 
Mapping/${sample}.filtered1.bam
done

#Viewing the output of bam file filtering
for sample in `cat raw_data/list.txt`
do
    samtools flagstat Mapping/${sample}.filtered1.bam
done

#Removing duplicates from a file: use the command rmdup
for sample in `cat raw_data/list.txt`
do 
   samtools rmdup Mapping/${sample}.sorted.bam  Mapping/${sample}.rdup 
done


#Left align BAM
for sample in `cat raw_data/list.txt`
do      
        cat Mapping/${sample}.rdup  | bamleftalign -f 
raw_data/hg19.chr5_12_17.fa -m 5 -c > Mapping/${sample}.leftAlign.bam

done
#-c - compressed, -m - max-iterations

#Recalibrating read mapping qualities
for sample in `cat raw_data/list.txt`
do
        samtools calmd -@ 32 -b Mapping/${sample}.leftAlign.bam 
raw_data/hg19.chr5_12_17.fa > Mapping/${sample}.recalibrate.bam
done

#Refiltering read mapping qualities
for sample in `cat raw_data/list.txt`
do
        bamtools filter -in Mapping/${sample}.recalibrate.bam -mapQuality 
"<=254" > Mapping/${sample}.refilter.bam
done

#installing VarScan
wget https://sourceforge.net/projects/varscan/files/VarScan.v2.3.9.jar	

#Converting data to pileup
mkdir Variants

for sample in `cat raw_data/list.txt`
do
        samtools mpileup -f raw_data/hg19.chr5_12_17.fa 
Mapping/${sample}.refilter.bam --min-MQ 1 --min-BQ 28 \
                > Variants/${sample}.pileup
done


#Call variants
java -jar VarScan.v2.3.9.jar somatic Variants/SLGFSK-N_231335.pileup \
        Variants/SLGFSK-T_231336.pileup Variants/SLGFSK \
        --normal-purity 1  --tumor-purity 0.5 --output-vcf 1 

#Merge vcf
bgzip Variants/SLGFSK.snp.vcf > Variants/SLGFSK.snp.vcf.gz
bgzip Variants/SLGFSK.indel.vcf > Variants/SLGFSK.indel.vcf.gz
tabix Variants/SLGFSK.snp.vcf.gz
tabix Variants/SLGFSK.indel.vcf.gz
bcftools merge Variants/SLGFSK.snp.vcf.gz Variants/SLGFSK.indel.vcf.gz 
--force-samples > Variants/SLGFSK.vcf

#Variant Annotation
#download jar file
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip

# Unzip file
unzip snpEff_latest_core.zip
		
#download snpEff database
java -jar  snpEff/snpEff.jar download hg19

#annotate variants
java -Xmx8g -jar snpEff/snpEff.jar hg19 Variants/SLGFSK.vcf > 
Variants/SLGFSK.ann.vcf

			

