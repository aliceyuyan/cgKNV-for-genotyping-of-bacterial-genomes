# cgKNV (core-genome based k-mer natural vector):
# A novel core genome clustering method proposed for population structure analysis and genotyping of bacterial genomes

This repository contains a complete, reproducible workflow for cgKNV by applied on 9 *B. anthracis* genomes (for example), as well as several process documents gotten from 346 *E. faecium*, 583 *B. anthracis*, and 1786 *M. abscessus* shown in the manuscript.  

The complete documents for the workflow of cgKNV applied on 9 *B. anthracis* genomes are shown in the directory and sub-directories of anthracis (~/anthracis), except several files that are too big to upload into the repository, in which the information of 9 *B. anthracis* genomes is shown in the file of ~/anthracis/Accession_number_list.xlsx.
![workflow of cgKNV](https://github.com/user-attachments/files/27337727/cgKNV_analysis_pipeline.pdf)

# The workflow of cgKNV contains the following steps:
⚙️ Step 1. Data acquisition using ascp or wget;

1.1 Single-end sequencing data
```bash
cd /home/alice/data/anthracis/test_single_end
for id in `cat Download_single_end_list.txt`; do ascp -P 33001 
  -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
  -QT -l 500m -k 1 -d "era-fasp@$id" ./ ; \
done
```
1.2 Paired-end sequencing data
```bash
cd /home/alice/data/anthracis/test_paired_end
for id in `cat Download_paired_end_list1.txt`; do \
  ascp -P 33001 \
  -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
  -QT -l 500m -k 1 -d "era-fasp@$id" ./ ; \
done
```
1.3 Assembled genome data
```bash
cd /home/alice/data/anthracis/test_genome_assembly
cat Download_assembly_list.txt | while read ID; do wget $ID; done
```

🌌 Step 2. Quality-control using Trimmomatic;

2.1 Single-end reads (.fastq.gz)
```bash
cd /home/alice/data/anthracis/test_single_end/
mkdir trimmomatic
cat single_list_trim.txt | while read ID; do \
  trimmomatic SE \
  -threads 10 \
  -phred33 \
  $ID".fastq.gz" \
  "./trimmomatic/"$ID"_single.fastq.gz" \
  ILLUMINACLIP:/home/alice/software/trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10 \
  LEADING:3 SLIDINGWINDOW:4:15 TRAILING:3 MINLEN:30; \
done
```
2.2 Paired-end reads (.fastq.gz)
```bash
cd /home/alice/data/anthracis/test_paired_end
mkdir trimmomatic
cat paired_list_trim.txt | while read ID; do \
  trimmomatic PE \
  -threads 10 \
  -phred33 \
  $ID"_1.fastq.gz" $ID"_2.fastq.gz" \
  "./trimmomatic/"$ID"_paired_1.fastq.gz" \
  "./trimmomatic/"$ID"_unpaired_1.fastq.gz" \
  "./trimmomatic/"$ID"_paired_2.fastq.gz" \
  "./trimmomatic/"$ID"_unpaired_2.fastq.gz" \
  ILLUMINACLIP:/home/alice/software/trimmomatic-0.39/adapters/NexteraPE-PE.fa:2:30:10 \
  LEADING:20 TRAILING:20 MINLEN:30; \
done
```

🧬 Step 3. SNP calling and consensus genomes generation using Snippy in the snippy-multi model;

3.1 Single-end-trimmed data

```bash
cd /home/alice/data/anthracis/test_single_end/trimmomatic/result
snippy-multi test_single.txt \
  --ref /home/alice/data/anthracis/NC_007530.2.fasta \
  --cpus 2 > test_single.sh
sh test_single.sh
```
3.2 Paired-end-trimmed data
```bash
cd /home/alice/data/anthracis/test_paired_end/trimmomatic/result
snippy-multi test_paired.txt \
  --ref /home/alice/data/anthracis/NC_007530.2.fasta \
  --cpus 2 > test_paired.sh
sh test_paired.sh
```
3.3 Assembled genomes
```bash
cd /home/alice/data/anthracis/test_assembled_genome
# Decompress
gunzip -d *_genomic.fna.gz
ls | grep _genomic.fna > list1
cat list1 | while read var; do echo ${var:0:15}; done > assembly_list.txt
paste assembly_list.txt list1 > list2
mkdir -p update
cat list2 | while read i j; do echo -e "$i\t/home/alice/data/anthracis/test_assembled_genome/$j" done > ./update/test_assembly.txt
cd /home/alice/data/anthracis/test_assembled_genome/update
snippy-multi test_assembly.txt \
  --ref /home/alice/data/anthracis/NC_007530.2.fasta \
  --cpus 2 > test_assembly.sh
sh test_assembly.sh
```
🪐 Step 4. Standardizing consensus genomes as in the .fasta file format;

4.1 Creating directory for consensus genomes
```bash
cd /home/alice/data/anthracis/
mkdir -p test_all_in_trimmomatic
```
4.2 Collecting single-end consensus genomes
```bash
cd /home/alice/data/anthracis/test_single_end/trimmomatic/result
for ID in $(cat ../../single_list_trim.txt); do
  cp "$ID/snps.consensus.subs.fa" "/home/alice/data/anthracis/test_all_in_trimmomatic/$ID.consensus.subs.fa"
done
```
4.3 Collecting paired-end consensus genomes
```bash
cd /home/alice/data/anthracis/test_paired_end/trimmomatic/result
for ID in $(cat ../../paired_list_trim.txt); do
  cp "$ID/snps.consensus.subs.fa" "/home/alice/data/anthracis/test_all_in_trimmomatic/$ID.consensus.subs.fa"
done
```
4.4 Collecting assembled genomes
```bash
cd /home/alice/data/anthracis/test_assembled_genome/update
for ID in $(cat ../assembly_list.txt); do
  cp "$ID/snps.consensus.subs.fa" "/home/alice/data/anthracis/test_all_in_trimmomatic/$ID.consensus.subs.fa"
done
```
4.5 Padding >Header to the consensus genomes
```bash
cd /home/alice/data/anthracis/test_all_in_trimmomatic
mkdir -p update
for i in $(cat Accession_ID_list.txt); do
  echo ">$i" > "update/$i.consensus.subs.fa"
  tail -n +2 "$i.consensus.subs.fa" >> "update/$i.consensus.subs.fa"
done
```
🏋️ Step 5. Bacterial genomes characterized using cgKNV and calculating distances of cgKNVs for bacterial genomes with the Hamming distance measure at ~/anthracis/test_all_in_trimmomatic/update;

Running cgKNV_analysis.m (for Windows)\
Running {the complete path of Matlab}/matlab  --nosplash   --nodesktop  cgKNV_analysis (for Linux)\
Output file: cgKNV_distance_matrix_hamming.meg.

🌲 Step 6. Clustering results visulization as NJ tree using MEGA/iTOL.

Input file: cgKNV_distance_matrix_hamming.meg.\
Output files: NJ tree (cgKNV_NJTree_hamming.pdf) and cgKNV_NJTree_hamming.nwk for visualized with iTOL (if necessary).
