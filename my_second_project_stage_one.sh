#downloading DNA.fa
wget https://raw.githubusercontent.com/HackBio-Internship/wale-home-tasks/main/DNA.fa

#Counting the number of sequences in DNA.fa
grep "^>" -c DNA.fa

# the total A, T, G & C counts for all the sequences in DNA.fa
awk -F "" 'BEGIN {totA=0; totT=0; totC=0; totG=0} !/^>/ {nA=gsub(/A/,A,$0); nT=gsub(/T/,T,$0); nC=gsub(/C/,C,$0); nG=gsub(/G/,G,$0); totA+=nA; totT+=nT; totC+=nC; totG+=nG} END {print "A-"totA"; T-"totT"; G-"totG"; C-"totC}' DNA.fa

#Setting up miniconda environment
sudo wget -c https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
sudo chmod +x Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x 86_64.sh
conda create --name HackBio_Workspace
conda activate HackBio_Workspace

# dowloading fastqc, fastp and seqtk
conda install -c bioconda fastqc  
conda install -c bioconda seqtk
conda install -c bioconda fastp

#dowloading datasets
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/ACBarrie_R1.fastq.gz?raw=true -O ACBarrie_R1.fastq.gz
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/ACBarrie_R2.fastq.gz?raw=true -O ACBarrie_R2.fastq.gz
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/Alsen_R1.fastq.gz?raw=true -O Alsen_R1.fastq.gz
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/Alsen_R2.fastq.gz?raw=true -O Alsen_R2.fastq.gz
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/Baxter_R1.fastq.gz?raw=true -O Baxter_R1.fastq.gz
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/Baxter_R2.fastq.gz?raw=true -O Baxter_R2.fastq.gz
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/Chara_R1.fastq.gz?raw=true -O Chara_R1.fastq.gz
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/Chara_R2.fastq.gz?raw=true -O Chara_R2.fastq.gz
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/Drysdale_R1.fastq.gz?raw=true -O Drysdale_R1.fastq.gz
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/Drysdale_R2.fastq.gz?raw=true -O Drysdale_R2.fastq.gz

mkdir output

#Perfoming fastqc
fastqc ACBarrie_R1.fastq.gz ACBarrie_R2.fastq.gz -output
fastqc Alsen_R1.fastq.gz Alsen_R2.fastq.gz -output
fastqc Baxter_R1.fastq.gz Baxter_R2.fastq.gz -output
fastqc Chara_R2.fastq.gz Chara.fastq.gz -output
fastqc Drysdale_R1.fastq.gz Drysdale_R2.fastq.gz -output

# Trimming adapters using fastp
fastp -i ACBarrie_R1.fastq.gz -I ACBarrie_R2.fastq.gz -o output/ACBarrie_R1_trimmed.fastq.gz -O output/ACBarrie_R2_trimmed.fastq.gz
fastp -i Alsen_R1.fastq.gz -I Alsen_R2.fastq.gz -o output/Alsen_R1_trimmed.fastq.gz -O output/Alsen_R2_trimmed.fastq.gz
fastp -i Baxter_R1.fastq.gz -I Baxter_R2.fastq.gz -o output/Baxter_R1_trimmed.fastq.gz -O output/Baxter_R2_trimmed.fastq.gz
fastp -i Chara_R1.fastq.gz -I Chara_R2.fastq.gz -o output/Chara_R1_trimmed.fastq.gz -O output/Chara_R2_trimmed.fastq.gz
fastp -i Drysdale_R1.fastq.gz -I Drysdale_R2.fastq.gz -o output/Drysdale_trimmed.fastq.gz -O output/Drysdale_R2_trimmed.fastq.gz

#using seqk to convert fastq files to fasta
seqtk seq -a ACBarrie_R1.fastq.gz > output/ACBarrie_R1.fa
seqtk seq -a ACBarrie_R2.fastq.gz > output/ACBarrie_R2.fa
seqtk seq -a Alsen_R1.fastq.gz > output/Alsen_R1.fa
seqtk seq -a Alsen_R2.fastq.gz > output/Alsen_R2.fa
seqtk seq -a Baxter_R1.fastq.gz > output/Baxter_R1.fa
seqtk seq -a Baxter_R2.fastq.gz > output/Baxter_R2.fa
seqtk seq -a Chara_R1.fastq.gz > output/Chara_R1.fa
seqtk seq -a Chara_R2.fastq.gz > output/Chara_R2.fa
seqtk seq -a Drysdale_R1.fastq.gz > output/Drysdale_R1.fa
seqtk seq -a Drysdale_R2.fastq.gz > output/Drysdale_R2.fa
