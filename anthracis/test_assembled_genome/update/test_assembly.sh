snippy --outdir 'GCA_000007845.1' --ctgs '/home/alice/data/anthracis/test_genome_assembly/GCA_000007845.1_ASM784v1_genomic.fna' --ref /home/alice/data/anthracis/NC_007530.2.fasta --cpus 12
snippy --outdir 'GCA_000008165.1' --ctgs '/home/alice/data/anthracis/test_genome_assembly/GCA_000008165.1_ASM816v1_genomic.fna' --ref /home/alice/data/anthracis/NC_007530.2.fasta --cpus 12
snippy --outdir 'GCF_001277955.1' --ctgs '/home/alice/data/anthracis/test_genome_assembly/GCF_001277955.1_ASM127795v1_genomic.fna' --ref /home/alice/data/anthracis/NC_007530.2.fasta --cpus 12
snippy-core --ref 'GCA_000007845.1/ref.fa' GCA_000007845.1 GCA_000008165.1 GCF_001277955.1
