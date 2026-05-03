# cgKNV (core-genome based k-mer natural vector):
# A novel core genome clustering method proposed for population structure analysis and genotyping of bacterial genomes.

# This repository contains a complete, reproducible workflow for cgKNV by applied on 9 B. anthracis genomes, as well as several process documents and the obtaining results from 346 E. faecium, 583 B. anthracis, and 1786 M. abscessus shown in the manuscript.  

The complete documents for the workflow of cgKNV applied on 9 B. anthracis genomes are shown in the directory of anthracis (~/anthracis), except several files that are too big to upload into the repository, which includes sequencing data of single_end (ERR1596541, ERR1596542, ERR1596543), paired_end (DRR014735, DRR014736, DRR014737), and assembled genomes (GCA_000007845.1, GCA_000008165.1, GCF_001277955.1), the informaiton of which are shown in the file of ~/anthracis/Accession_number_list.xlsx.

The workflow for B. anthracis is composed of following steps (Steps 1-6):

---
## 📂 Project Structure
The directories of cgKNV are shown as follows:
~/anthracis/
├── NC_007530.2.fasta              # Reference genome
├── Core_gene_positions.txt        # Core gene coordinates annotated in the reference genome NC_007530.2.fasta
├── test_single_end/               # Single-end raw data
├── test_paired_end/               # Paired-end raw data
├── test_assembled_genome/         # Assembled genomes
└── test_all_in_trimmomatic/       # Consolidated consensus sequences

⚙️ 1. Raw Data Acquisition
Prerequisites
Tools: Aspera Connect (ascp), wget
Input Lists: Download_single_end_list.txt, Download_paired_end_list1.txt, etc.
Commands (Unmodified)
Single-end Sequencing:

```bash
cd /home/alice/data/anthracis/test_single_end
for id in `cat Download_single_end_list.txt`; do \
  ascp -P 33001 \
  -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
  -QT -l 500m -k 1 -d "era-fasp@$id" ./ ; \
done

Paired-end Sequencing:

```bash
cd /home/alice/data/anthracis/test_paired_end
for id in `cat Download_paired_end_list1.txt`; do \
  ascp -P 33001 \
  -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
  -QT -l 500m -k 1 -d "era-fasp@$id" ./ ; \
done

Genome Assembly (NCBI):

```bash
cd /home/alice/data/anthracis/test_genome_assembly
cat Download_assembly_list.txt | while read ID; do wget $ID; done

🧹 2. Quality Control with Trimmomatic
Single-end Reads

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

Paired-end Reads

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

🧬 3. SNP Calling & Standardization
Single-end Data

```bash
cd /home/alice/data/anthracis/test_single_end/trimmomatic/result
snippy-multi test_single.txt \
  --ref /home/alice/data/anthracis/NC_007530.2.fasta \
  --cpus 2 > test_single.sh
sh test_single.sh

Paired-end Data

```bash
cd /home/alice/data/anthracis/test_paired_end/trimmomatic/result
snippy-multi test_paired.txt \
  --ref /home/alice/data/anthracis/NC_007530.2.fasta \
  --cpus 2 > test_paired.sh
sh test_paired.sh

Assembled Genomes

```bash
cd /home/alice/data/anthracis/test_assembled_genome
# Decompress
gunzip -d *_genomic.fna.gz
# Generate list (Note: Adjust substring length if needed)
ls | grep _genomic.fna > list1
cat list1 | while read var; do echo ${var:0:15}; done > assembly_list.txt
paste assembly_list.txt list1 > list2
mkdir -p update
cat list2 | while read i j; do echo -e "$i\t/home/alice/data/anthracis/test_assembled_genome/$j" done > ./update/test_assembly.txt

# Run Snippy
cd /home/alice/data/anthracis/test_assembled_genome/update
snippy-multi test_assembly.txt \
  --ref /home/alice/data/anthracis/NC_007530.2.fasta \
  --cpus 2 > test_assembly.sh
sh test_assembly.sh

Standardize Consensus Sequences

```bash
# 1. Create Consolidated Directory
cd /home/alice/data/anthracis/
mkdir -p test_all_in_trimmomatic

# 2. Collect Files (Single-end)
cd /home/alice/data/anthracis/test_single_end/trimmomatic/result
for ID in $(cat ../../single_list_trim.txt); do
  cp "$ID/snps.consensus.subs.fa" "/home/alice/data/anthracis/test_all_in_trimmomatic/$ID.consensus.subs.fa"
done

# 3. Collect Files (Paired-end)
cd /home/alice/data/anthracis/test_paired_end/trimmomatic/result
for ID in $(cat ../../paired_list_trim.txt); do
  cp "$ID/snps.consensus.subs.fa" "/home/alice/data/anthracis/test_all_in_trimmomatic/$ID.consensus.subs.fa"
done

# 4. Collect Files (Assembled)
cd /home/alice/data/anthracis/test_assembled_genome/update
for ID in $(cat ../assembly_list.txt); do
  cp "$ID/snps.consensus.subs.fa" "/home/alice/data/anthracis/test_all_in_trimmomatic/$ID.consensus.subs.fa"
done

# 5. Standardize Headers
cd /home/alice/data/anthracis/test_all_in_trimmomatic
mkdir -p update
for i in $(cat Accession_ID_list.txt); do
  echo ">$i" > "update/$i.consensus.subs.fa"
  tail -n +2 "$i.consensus.subs.fa" >> "update/$i.consensus.subs.fa"
done

🧮 4. cgKNV Vector Calculation (MATLAB)
Tool: MATLAB (R2025b)
Copy the standardized consensus files (update/*.consensus.subs.fa) and Core_gene_positions.txt to MATLAB.
Run the script to extract core genes, determine k-value, and compute KNVs.
Output: cgKNV_distance_matrix_hamming.meg

🌲 5. Phylogenetic Analysis
MEGA11
Import cgKNV_distance_matrix_hamming.meg.
Construct a Neighbor-Joining (NJ) Tree.
Export as cgKNV_NJTree_hamming.nwk.
iTOL
Upload the cgKNV_NJTree_hamming.nwk file.
Customize and export as cgKNV_NJTree_hamming.pdf.
