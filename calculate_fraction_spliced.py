import pysam
from collections import Counter
from CGAT import GTF, IOTools
import numpy
import sys

# import files 
print(sys.argv)
bamfile = pysam.AlignmentFile(sys.argv[1], "rb")
gtffile=sys.argv[2]
outFile = open(sys.argv[3], 'w')

outFile.write(str('GeneID')+'\t'+str('TxID')+'\t'+str('Strand')+'\t'+str('Chr')+'\t'+str('intron_coor')+'\t'+str('intron_length')+'\t'+str('intron_num')+'\t'+str('Fraction_3p_spliced')+'\t'+str('Fraction_3p_annotated_splice')+'\t'+str('Exon_Exon_3p_counts')+'\t'+str('Unspliced_3p_counts')+'\t'+str('spliced_3p_uncounted')+'\t'+str('spliced_3p_block_too_small')+'\t'+str('unspliced_3p_uncounted')+'\t'+str('antisense_3p_counts')+'\t'+str('Fraction_5p_spliced')+'\t'+str('Fraction_5p_annotated_splice')+'\t'+str('Exon_Exon_5p_counts')+'\t'+str('Unspliced_5p_counts')+'\t'+str('spliced_5p_uncounted')+'\t'+str('spliced_5p_block_too_small')+'\t'+str('unspliced_5p_uncounted')+'\t'+str('antisense_5p_counts')+'\n')

for transcript in GTF.transcript_iterator(GTF.iterator(IOTools.openFile(gtffile))):
   introns = GTF.toIntronIntervals(transcript)
   txid=transcript[0].transcript_id
   geneid=transcript[0].gene_id
   strand=transcript[0].strand
   chr=transcript[0].contig
   min_overlap=10
   searchwin=0
   if strand == "+":
      intron_num=0
   elif strand == "-":
      intron_num=len(introns)+1
   for intron in introns:
      if strand == "+":
         intron5p=intron[0]
         intron3p=intron[1]
      elif strand == "-":
         intron5p=intron[1]
         intron3p=intron[0]
      intron_length=intron[1]-intron[0]
      if strand == "+":
         intron_num += 1
      elif strand == "-":
         intron_num += -1
      Exon_Exon_3p=0
      spliced_3p_uncounted=0
      spliced_3p_too_small_block=0
      unspliced3p=0
      unspliced3p_uncounted=0
      antisense_3p_counts=0
      Exon_Exon_5p=0
      spliced_5p_uncounted=0
      spliced_5p_too_small_block=0
      unspliced5p=0
      unspliced5p_uncounted=0
      antisense_5p_counts=0
      fraction_3p_spliced=0
      frac_3p_anno_spliced=0
      fraction_5p_spliced=0
      frac_5p_anno_spliced=0
      found = set()
      reads = bamfile.fetch(contig=chr, start=intron3p-1-searchwin, stop=intron3p+searchwin)
      for read in reads:
         #select sense reads
         if (strand == "+" and read.is_read1 and read.is_reverse) or (strand == "+" and read.is_read2 and not read.is_reverse) or (strand == "-" and read.is_read1 and not read.is_reverse) or (strand == "-" and read.is_read2 and read.is_reverse):
            #get reads with N in cigar string
            if 'N' in read.cigarstring:
               blocks = read.get_blocks()
               starts, ends = zip(*blocks)
               # chuck reads with tiny blocks, need all blocks to be >= 10
               if numpy.amin(numpy.subtract((ends), (starts))) >= min_overlap:
                  if intron[0] in ends and intron[1] in starts and read.reference_end-intron[1]>=min_overlap and intron[0]-read.reference_start>=min_overlap:
                     readname = str(read.query_name)+"_3SS_"+str(intron3p)
                     if readname not in found:
                        found.add(readname)
                        Exon_Exon_3p += 1
                  else:
                     readname = str(read.query_name)+"_3SS_"+str(intron3p)
                     if readname not in found:
                        found.add(readname)
                        spliced_3p_uncounted += 1
               else:
                  readname = str(read.query_name)+"_3SS_"+str(intron3p)
                  if readname not in found:
                     spliced_3p_too_small_block += 1
            elif (read.reference_start <= intron3p - min_overlap and read.reference_end >= intron3p + min_overlap):
               readname = str(read.query_name)+"_3SS_"+str(intron3p)
               if readname not in found:
                  found.add(readname)
                  unspliced3p += 1
            else:
               readname = str(read.query_name)+"_3SS_"+str(intron3p)
               if readname not in found:
                  unspliced3p_uncounted += 1
         else:
            readname = str(read.query_name)+"_3SS_"+str(intron3p)
            if readname not in found:
               found.add(readname)
               antisense_3p_counts += 1
         if (Exon_Exon_3p + unspliced3p) > 0:
            fraction_3p_spliced= Exon_Exon_3p / (Exon_Exon_3p + unspliced3p)
         else:
            fraction_3p_spliced="NA"
         if (Exon_Exon_3p + spliced_3p_uncounted) > 0:
            frac_3p_anno_spliced= Exon_Exon_3p / (Exon_Exon_3p + spliced_3p_uncounted)
         else:
            frac_3p_anno_spliced="NA"
      
      ######## Repeat for the 5' SS
      reads5p = bamfile.fetch(contig=chr, start=intron5p-1-searchwin, stop=intron5p+searchwin)
      for read in reads5p:
         #select sense reads
         if (strand == "+" and read.is_read1 and read.is_reverse) or (strand == "+" and read.is_read2 and not read.is_reverse) or (strand == "-" and read.is_read1 and not read.is_reverse) or (strand == "-" and read.is_read2 and read.is_reverse):
            #get reads with N in cigar string
            if 'N' in read.cigarstring:
               blocks = read.get_blocks()
               starts, ends = zip(*blocks)
               # chuck reads with tiny blocks, need all blocks to be >= 10
               if numpy.amin(numpy.subtract((ends), (starts))) >= min_overlap:
                  if intron[0] in ends and intron[1] in starts and read.reference_end-intron[1]>=min_overlap and intron[0]-read.reference_start>=min_overlap:
                     readname = str(read.query_name)+"_5SS_"+str(intron5p)
                     if readname not in found:
                        found.add(readname)
                        Exon_Exon_5p += 1
                  else:
                     readname = str(read.query_name)+"_5SS_"+str(intron5p)
                     if readname not in found:
                        found.add(readname)
                        spliced_5p_uncounted += 1
               else:
                  readname = str(read.query_name)+"_5SS_"+str(intron5p)
                  if readname not in found:
                     spliced_5p_too_small_block += 1
            elif (read.reference_start <= intron5p - min_overlap and read.reference_end >= intron5p + min_overlap):
               readname = str(read.query_name)+"_5SS_"+str(intron5p)
               if readname not in found:
                  found.add(readname)
                  unspliced5p += 1
            else:
               readname = str(read.query_name)+"_5SS_"+str(intron5p)
               if readname not in found:
                  unspliced5p_uncounted += 1
         else:
            readname = str(read.query_name)+"_5SS_"+str(intron5p)
            if readname not in found:
               found.add(readname)
               antisense_5p_counts += 1
         if (Exon_Exon_5p + unspliced5p) > 0:
            fraction_5p_spliced= Exon_Exon_5p / (Exon_Exon_5p + unspliced5p)
         else:
            fraction_5p_spliced="NA"
         if (Exon_Exon_5p + spliced_5p_uncounted) > 0:
            frac_5p_anno_spliced= Exon_Exon_5p / (Exon_Exon_5p + spliced_5p_uncounted)
         else:
            frac_5p_anno_spliced="NA"
      ### output to file
      outFile.write(str(geneid)+'\t'+str(txid)+'\t'+str(strand)+'\t'+str(chr)+'\t'+str(intron)+'\t'+str(int(intron_length))+'\t'+str(int(intron_num))+'\t'+str(fraction_3p_spliced)+'\t'+str(frac_3p_anno_spliced)+'\t'+str(Exon_Exon_3p)+'\t'+str(unspliced3p)+'\t'+str(spliced_3p_uncounted)+'\t'+str(spliced_3p_too_small_block)+'\t'+str(unspliced3p_uncounted)+'\t'+str(antisense_3p_counts)+'\t'+str(fraction_5p_spliced)+'\t'+str(frac_5p_anno_spliced)+'\t'+str(Exon_Exon_5p)+'\t'+str(unspliced5p)+'\t'+str(spliced_5p_uncounted)+'\t'+str(spliced_5p_too_small_block)+'\t'+str(unspliced5p_uncounted)+'\t'+str(antisense_5p_counts)+'\n')


# close all files
bamfile.close()
outFile.close()