#############################################
##### Sequencing the mutations analysis #####
#############################################

# Mapping to the reference genome
minimap2 -ax map-ont -Y GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna raw_reads/T624WQ_1_NF_12.fastq > NF_12.sam
minimap2 -ax map-ont -Y GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna raw_reads/T624WQ_2_NF_46.fastq > NF_46.sam
minimap2 -ax map-ont -Y GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna raw_reads/T624WQ_3_NF_48.fastq > NF_48.sam
minimap2 -ax map-ont -Y GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna raw_reads/T624WQ_4_NF_58.fastq > NF_58.sam
minimap2 -ax map-ont -Y GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna raw_reads/T624WQ_5_NF_61.fastq > NF_61.sam
minimap2 -ax map-ont -Y GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna raw_reads/T624WQ_6_NF_64.fastq > NF_64.sam
minimap2 -ax map-ont -Y GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna raw_reads/T624WQ_7_NF_69.fastq > NF_69.sam
minimap2 -ax map-ont -Y GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna raw_reads/T624WQ_8_NF_72.fastq > NF_72.sam
minimap2 -ax map-ont -Y GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna raw_reads/T624WQ_9_NF_76.fastq > NF_76.sam
minimap2 -ax map-ont -Y GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna raw_reads/T624WQ_10_NF_78.fastq > NF_78.sam
minimap2 -ax map-ont -Y GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna raw_reads/T624WQ_11_NF_80.fastq > NF_80.sam

# Converting sam to bam
samtools view -@ 8 -b -o NF_12.bam NF_12.sam
samtools view -@ 8 -b -o NF_46.bam NF_46.sam
samtools view -@ 8 -b -o NF_48.bam NF_48.sam
samtools view -@ 8 -b -o NF_58.bam NF_58.sam
samtools view -@ 8 -b -o NF_61.bam NF_61.sam
samtools view -@ 8 -b -o NF_64.bam NF_64.sam
samtools view -@ 8 -b -o NF_69.bam NF_69.sam
samtools view -@ 8 -b -o NF_72.bam NF_72.sam
samtools view -@ 8 -b -o NF_76.bam NF_76.sam
samtools view -@ 8 -b -o NF_78.bam NF_78.sam
samtools view -@ 8 -b -o NF_80.bam NF_80.sam

# Sorting
samtools sort NF_12.bam -o NF_12.sorted.bam
samtools sort NF_46.bam -o NF_46.sorted.bam
samtools sort NF_48.bam -o NF_48.sorted.bam
samtools sort NF_58.bam -o NF_58.sorted.bam
samtools sort NF_61.bam -o NF_61.sorted.bam
samtools sort NF_64.bam -o NF_64.sorted.bam
samtools sort NF_69.bam -o NF_69.sorted.bam
samtools sort NF_72.bam -o NF_72.sorted.bam
samtools sort NF_76.bam -o NF_76.sorted.bam
samtools sort NF_78.bam -o NF_78.sorted.bam
samtools sort NF_80.bam -o NF_80.sorted.bam

# Filtering reads based on region
samtools view -b -o NF_12.sorted.filtered.bam NF_12.sorted.bam NT_033777.3:5465755-5466215
samtools view -b -o NF_46.sorted.filtered.bam NF_46.sorted.bam NT_033777.3:5465755-5466215
samtools view -b -o NF_48.sorted.filtered.bam NF_48.sorted.bam NT_033777.3:5465755-5466215
samtools view -b -o NF_58.sorted.filtered.bam NF_58.sorted.bam NT_033777.3:5465755-5466215
samtools view -b -o NF_61.sorted.filtered.bam NF_61.sorted.bam NT_033777.3:5465755-5466215
samtools view -b -o NF_64.sorted.filtered.bam NF_64.sorted.bam NT_033777.3:5465755-5466215
samtools view -b -o NF_69.sorted.filtered.bam NF_69.sorted.bam NT_033777.3:5465755-5466215
samtools view -b -o NF_72.sorted.filtered.bam NF_72.sorted.bam NT_033777.3:5465755-5466215
samtools view -b -o NF_76.sorted.filtered.bam NF_76.sorted.bam NT_033777.3:5465755-5466215
samtools view -b -o NF_78.sorted.filtered.bam NF_78.sorted.bam NT_033777.3:5465755-5466215
samtools view -b -o NF_80.sorted.filtered.bam NF_80.sorted.bam NT_033777.3:5465755-5466215

# Indexing for IGV and whatshap
samtools faidx GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna
samtools index NF_12.sorted.filtered.bam
samtools index NF_46.sorted.filtered.bam
samtools index NF_48.sorted.filtered.bam
samtools index NF_58.sorted.filtered.bam
samtools index NF_61.sorted.filtered.bam
samtools index NF_64.sorted.filtered.bam
samtools index NF_69.sorted.filtered.bam
samtools index NF_72.sorted.filtered.bam
samtools index NF_76.sorted.filtered.bam
samtools index NF_78.sorted.filtered.bam
samtools index NF_80.sorted.filtered.bam

# Variant calling
freebayes -f GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna NF_12.sorted.filtered.bam >NF_12.vcf
freebayes -f GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna NF_46.sorted.filtered.bam >NF_46.vcf
freebayes -f GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna NF_48.sorted.filtered.bam >NF_48.vcf
freebayes -f GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna NF_58.sorted.filtered.bam >NF_58.vcf
freebayes -f GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna NF_61.sorted.filtered.bam >NF_61.vcf
freebayes -f GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna NF_64.sorted.filtered.bam >NF_64.vcf
freebayes -f GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna NF_69.sorted.filtered.bam >NF_69.vcf
freebayes -f GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna NF_72.sorted.filtered.bam >NF_72.vcf
freebayes -f GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna NF_76.sorted.filtered.bam >NF_76.vcf
freebayes -f GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna NF_78.sorted.filtered.bam >NF_78.vcf
freebayes -f GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna NF_80.sorted.filtered.bam >NF_80.vcf

# Phasing haplotypes
whatshap phase -o NF_12.phased.vcf --reference=GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna NF_12.vcf NF_12.sorted.filtered.bam
whatshap phase -o NF_46.phased.vcf --reference=GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna NF_46.vcf NF_46.sorted.filtered.bam
whatshap phase -o NF_48.phased.vcf --reference=GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna NF_48.vcf NF_48.sorted.filtered.bam
whatshap phase -o NF_58.phased.vcf --reference=GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna NF_58.vcf NF_58.sorted.filtered.bam
whatshap phase -o NF_61.phased.vcf --reference=GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna NF_61.vcf NF_61.sorted.filtered.bam
whatshap phase -o NF_64.phased.vcf --reference=GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna NF_64.vcf NF_64.sorted.filtered.bam
whatshap phase -o NF_69.phased.vcf --reference=GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna NF_69.vcf NF_69.sorted.filtered.bam
whatshap phase -o NF_72.phased.vcf --reference=GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna NF_72.vcf NF_72.sorted.filtered.bam
whatshap phase -o NF_76.phased.vcf --reference=GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna NF_76.vcf NF_76.sorted.filtered.bam
whatshap phase -o NF_78.phased.vcf --reference=GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna NF_78.vcf NF_78.sorted.filtered.bam
whatshap phase -o NF_80.phased.vcf --reference=GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna NF_80.vcf NF_80.sorted.filtered.bam

# Compress vcf files
bgzip NF_12.phased.vcf
bgzip NF_46.phased.vcf
bgzip NF_48.phased.vcf
bgzip NF_58.phased.vcf
bgzip NF_61.phased.vcf
bgzip NF_64.phased.vcf
bgzip NF_69.phased.vcf
bgzip NF_72.phased.vcf
bgzip NF_76.phased.vcf
bgzip NF_78.phased.vcf
bgzip NF_80.phased.vcf

# Indexing vcf files
tabix -p vcf NF_12.phased.vcf.gz
tabix -p vcf NF_46.phased.vcf.gz
tabix -p vcf NF_48.phased.vcf.gz
tabix -p vcf NF_58.phased.vcf.gz
tabix -p vcf NF_61.phased.vcf.gz
tabix -p vcf NF_64.phased.vcf.gz
tabix -p vcf NF_69.phased.vcf.gz
tabix -p vcf NF_72.phased.vcf.gz
tabix -p vcf NF_76.phased.vcf.gz
tabix -p vcf NF_78.phased.vcf.gz
tabix -p vcf NF_80.phased.vcf.gz



# Creating consensus haplotypes
bcftools consensus -H 1 -f GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna NF_12.phased.vcf.gz > NF_12_haplotype1.fasta
bcftools consensus -H 2 -f GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna NF_12.phased.vcf.gz > NF_12_haplotype2.fasta

bcftools consensus -H 1 -f GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna NF_46.phased.vcf.gz > NF_46_haplotype1.fasta
bcftools consensus -H 2 -f GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna NF_46.phased.vcf.gz > NF_46_haplotype2.fasta

bcftools consensus -H 1 -f GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna NF_48.phased.vcf.gz > NF_48_haplotype1.fasta
bcftools consensus -H 2 -f GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna NF_48.phased.vcf.gz > NF_48_haplotype2.fasta

bcftools consensus -H 1 -f GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna NF_58.phased.vcf.gz > NF_58_haplotype1.fasta
bcftools consensus -H 2 -f GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna NF_58.phased.vcf.gz > NF_58_haplotype2.fasta

bcftools consensus -H 1 -f GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna NF_61.phased.vcf.gz > NF_61_haplotype1.fasta
bcftools consensus -H 2 -f GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna NF_61.phased.vcf.gz > NF_61_haplotype2.fasta

bcftools consensus -H 1 -f GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna NF_64.phased.vcf.gz > NF_64_haplotype1.fasta
bcftools consensus -H 2 -f GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna NF_64.phased.vcf.gz > NF_64_haplotype2.fasta

bcftools consensus -H 1 -f GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna NF_69.phased.vcf.gz > NF_69_haplotype1.fasta
bcftools consensus -H 2 -f GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna NF_69.phased.vcf.gz > NF_69_haplotype2.fasta

bcftools consensus -H 1 -f GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna NF_72.phased.vcf.gz > NF_72_haplotype1.fasta
bcftools consensus -H 2 -f GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna NF_72.phased.vcf.gz > NF_72_haplotype2.fasta

bcftools consensus -H 1 -f GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna NF_76.phased.vcf.gz > NF_76_haplotype1.fasta
bcftools consensus -H 2 -f GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna NF_76.phased.vcf.gz > NF_76_haplotype2.fasta

bcftools consensus -H 1 -f GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna NF_78.phased.vcf.gz > NF_78_haplotype1.fasta
bcftools consensus -H 2 -f GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna NF_78.phased.vcf.gz > NF_78_haplotype2.fasta

bcftools consensus -H 1 -f GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna NF_80.phased.vcf.gz > NF_80_haplotype1.fasta
bcftools consensus -H 2 -f GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna NF_80.phased.vcf.gz > NF_80_haplotype2.fasta

# Selecting only the RpL35A region in the genome
samtools faidx GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna NT_033777.3:5465755-5466215 > Reference_RpL35A.fasta

samtools faidx NF_12_haplotype1.fasta NT_033777.3:5465755-5466215 > NF_12_haplotype1_RpL35A.fasta
samtools faidx NF_12_haplotype2.fasta NT_033777.3:5465755-5466215 > NF_12_haplotype2_RpL35A.fasta

samtools faidx NF_46_haplotype1.fasta NT_033777.3:5465755-5466215 > NF_46_haplotype1_RpL35A.fasta
samtools faidx NF_46_haplotype2.fasta NT_033777.3:5465755-5466215 > NF_46_haplotype2_RpL35A.fasta

samtools faidx NF_48_haplotype1.fasta NT_033777.3:5465755-5466215 > NF_48_haplotype1_RpL35A.fasta
samtools faidx NF_48_haplotype2.fasta NT_033777.3:5465755-5466215 > NF_48_haplotype2_RpL35A.fasta

samtools faidx NF_58_haplotype1.fasta NT_033777.3:5465755-5466215 > NF_58_haplotype1_RpL35A.fasta
samtools faidx NF_58_haplotype2.fasta NT_033777.3:5465755-5466215 > NF_58_haplotype2_RpL35A.fasta

samtools faidx NF_61_haplotype1.fasta NT_033777.3:5465755-5466215 > NF_61_haplotype1_RpL35A.fasta
samtools faidx NF_61_haplotype2.fasta NT_033777.3:5465755-5466215 > NF_61_haplotype2_RpL35A.fasta

samtools faidx NF_64_haplotype1.fasta NT_033777.3:5465755-5466215 > NF_64_haplotype1_RpL35A.fasta
samtools faidx NF_64_haplotype2.fasta NT_033777.3:5465755-5466215 > NF_64_haplotype2_RpL35A.fasta

samtools faidx NF_69_haplotype1.fasta NT_033777.3:5465755-5466215 > NF_69_haplotype1_RpL35A.fasta
samtools faidx NF_69_haplotype2.fasta NT_033777.3:5465755-5466215 > NF_69_haplotype2_RpL35A.fasta

samtools faidx NF_72_haplotype1.fasta NT_033777.3:5465755-5466215 > NF_72_haplotype1_RpL35A.fasta
samtools faidx NF_72_haplotype2.fasta NT_033777.3:5465755-5466215 > NF_72_haplotype2_RpL35A.fasta

samtools faidx NF_76_haplotype1.fasta NT_033777.3:5465755-5466215 > NF_76_haplotype1_RpL35A.fasta
samtools faidx NF_76_haplotype2.fasta NT_033777.3:5465755-5466215 > NF_76_haplotype2_RpL35A.fasta

samtools faidx NF_78_haplotype1.fasta NT_033777.3:5465755-5466215 > NF_78_haplotype1_RpL35A.fasta
samtools faidx NF_78_haplotype2.fasta NT_033777.3:5465755-5466215 > NF_78_haplotype2_RpL35A.fasta

samtools faidx NF_80_haplotype1.fasta NT_033777.3:5465755-5466215 > NF_80_haplotype1_RpL35A.fasta
samtools faidx NF_80_haplotype2.fasta NT_033777.3:5465755-5466215 > NF_80_haplotype2_RpL35A.fasta

# Printing the haplotypes to put them in a fasta file manually
cat NF_12_haplotype1_RpL35A.fasta
cat NF_12_haplotype2_RpL35A.fasta
cat NF_46_haplotype1_RpL35A.fasta
cat NF_46_haplotype2_RpL35A.fasta
cat NF_48_haplotype1_RpL35A.fasta
cat NF_48_haplotype2_RpL35A.fasta
cat NF_58_haplotype1_RpL35A.fasta
cat NF_58_haplotype2_RpL35A.fasta
cat NF_61_haplotype1_RpL35A.fasta
cat NF_61_haplotype2_RpL35A.fasta
cat NF_64_haplotype1_RpL35A.fasta
cat NF_64_haplotype2_RpL35A.fasta
cat NF_69_haplotype1_RpL35A.fasta
cat NF_69_haplotype2_RpL35A.fasta
cat NF_72_haplotype1_RpL35A.fasta
cat NF_72_haplotype2_RpL35A.fasta
cat NF_76_haplotype1_RpL35A.fasta
cat NF_76_haplotype2_RpL35A.fasta
cat NF_78_haplotype1_RpL35A.fasta
cat NF_78_haplotype2_RpL35A.fasta
cat NF_80_haplotype1_RpL35A.fasta
cat NF_80_haplotype2_RpL35A.fasta

# Manually fixing two phasing errors (after inspection with IGV)
NF_72_haplotype1
GGGG -> TGGC
NF_72_haplotype2
TGGC -> GGGG

NF_78_haplotype1
GGGG -> TGGC
NF_78_haplotype2
TGGC -> GGGG

# Reverse complemented with webtool: https://www.bioinformatics.org/sms2/rev_comp.html



