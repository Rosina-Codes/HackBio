first_name="Rosina"
last_name="Carr"
# version where all variables are printed on the same line
echo $first_name $last_name
# version where all variables are printed on the different lines
echo -e "$first_name\n$last_name"
#Creating a folder with my name
mkdir Carr
#navigating into it
cd Carr
#Creating a folder in Carr
mkdir biocomputing && cd biocomputing
#downloading files into the biocomuputing folder
wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.fna .
wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.gbk .
wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.gbk .
# moving downloaded files into the parent directory Carr
mv wildtype.fna wildtype.gbk wildtype.gbk.1 ..
# navigating into the parent directory
cd ..
# deleting the wildtype.gbk duplicate
rm wildtype.gbk.1
# double-checking the deletion of the duplicate file
ls
# finding the pattern "tatatata"
grep "tatatata" wildtype.fna
#printing al the line wit "tatatata" into a file
grep -n  "tatatata" wildtype.fna > mutant.fna
#clearing my terminal
clear
#listingthe files in the two folders 
ls -R

#Assignment 2
#installing figlet
sudo apt instal figlet
#displaying a graphical representation of my name using figlet
figlet 'Carr Rosina'
# creating and navigating into a folder called compare
mkdir compare && cd compare
#dowloading a file from the internet into my working directory
wget https://www.bioinformatics.babraham.ac.uk/training/Introduction%20to%20Unix/unix_intro_data.tar.gz
#unzipping the dowloaded file
gunzip unix_intro_data.tar.gz
# untarring the unzipped folder
tar -xvf unix_intro_data.tar
#navigating into seqmonk_genomes/Saccharomyces cerevisiae/EF4 
cd seqmonk_genomes/"Saccharomyces cerevisiae"/EF4
# finding rRNAs in the Mito.dat file
grep "rRNA" Mito.dat
#Copying Mito.dat to the compare directory
cp Mito.dat ../../..
# opening the Mito.dat using nano for editing
nano Mito.dat
# renaming the Mito.dat file
mv Mito.dat Mitochondrion.txt
#navigating into the FastQ_Data
cd FastQ_Data
#total number of lines in lane8_DD_P4_TTAGGC_L008_R1.fastq.gz
wc -l lane8_DD_P4_TTAGGC_L008_R1.fastq.gz
#total number of lines in all fastq.gz files and save it as a new file
wc -l *.fastq.gz > total_fastq.gz_line_number