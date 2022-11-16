# calculate_fraction_spliced

Python script to calculate the fraction of annotated spliced reads at intron-exon junctions in TT-seq, written by Benjamin J. E. Martin (Karen Adelman Lab). 

Inputs:
  1. BAM file (paired end reverse stranded BAM files only)
  2. GTF file for gene annotations 
  3. Output file name

## Running the script

Step 1: Load the necessary conda environment
```
module load gcc/6.2.0
module load python/3.6.0
module load conda2/4.2.13

## This script requires cgat-tools, install cgat-tools if needed
cgat-tools: https://github.com/CGATOxford/cgat/blob/master/install-CGAT-tools.sh

```

Step 2: Sort the GTF by transcript id (note, this is assuming the ENSEMBL ordering of the 9th column):

```
sort -t ";" -k 3 your_gtf.gtf > your_gtf_sorted_by_tx_id.gtf
```

Step 3: Run the script

python calculate_fraction_spliced.py {bam} {gtf} {output_filename}

```
python calculate_fraction_spliced.py your_data.bam your_gtf_sorted_by_tx_id.gtf output_file_name.txt

```

## The output file contains the following columns:

1. GeneID - Gene ID 
2. TxID - Transcript ID 
3. Strand - Strand
4. Chr - Chromosome
5. intron_coor - the intron coordinates
6. intron_length - the intron length
7. intron_num - the intron number (starting from the TSS)
8. Fraction_3p_spliced - fraction of annotated 3' spliced reads (vs spliced+unspliced): Exon_Exon_3p_counts / (Exon_Exon_3p_counts + Unspliced_3p_counts)
9. Fraction_3p_annotated_splice - the fraction of 3' spliced reads that are the annotated exon-exon splicing events: Exon_Exon_3p_counts / (Exon_Exon_3p_counts + spliced_3p_uncounted)
10. Exon_Exon_3p_counts - the number of reads counted across the annotated exon-exon splice site
11. Unspliced_3p_counts - the number of unspliced reads counted across the 3' splice site
12. spliced_3p_uncounted - the number of spliced reads, aligned to 3' exon, but not consisting of annotated exon-exon splicing event (either a different 5' or 3' splice site gets counted here)
13. spliced_3p_block_too_small - reads with detected splicing but less than 10bp aligned on both ends of splicing (either at the 3' splice site or at the other end of splice junction
14. unspliced_3p_uncounted - unspliced reads not counted because the read did not contain 10bp aligned on either end of 3' splice site
15. antisense_3p_counts - the number of 3' reads counted on the antisense strand
16. Fraction_5p_spliced- fraction of annotated 5' aligning reads that are spliced: Exon_Exon_5p_counts / (Exon_Exon_5p_counts + Unspliced_5p_counts)
17. Fraction_5p_annotated_splice - the fraction of 5' spliced reads that are the annotated exon-exon spliced: Exon_Exon_5p_counts / (Exon_Exon_5p_counts + spliced_5p_uncounted)
18. Exon_Exon_5p_counts - the number of reads counted across the annotated exon-exon splice site
19. Unspliced_5p_counts - the number of unspliced reads counted across the 5' splice site
20. spliced_5p_uncounted - the number of spliced reads, aligned to 5' exon, but not consisting of annotated exon-exon splicing
21. spliced_5p_block_too_small  - reads with detected splicing but less than 10bp aligned on both ends of splicing (either at the 5' splice site or at the other end of splice junction
22. unspliced_5p_uncounted - unspliced reads not counted because the read did not contain 10bp aligned on either end of 5' splice site
23. antisense_5p_counts - the number of 5' reads counted on the antisense strand

Additional Notes: 
1. Written for analysis of constitutive exons/introns
2. This script counts TT-seq fragments i.e., if both read1 and read2 align to the 3' splice site this will count as 1 read in the count totals.
