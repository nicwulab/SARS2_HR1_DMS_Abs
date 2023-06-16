#!/usr/bin/python
import os
import sys
import operator
from Bio import SeqIO
from collections import Counter

def ProcessMultilib(Rfile):
  print ("Reading %s" % Rfile)
  records = SeqIO.parse(Rfile,"fastq")
  variants = [] 
  record_count = 0
  for record in records:
    record_count += 1
    Rseq  = record.seq
    Rroi = Rseq
    if ((Rroi[0:24] == "ACATCTGCCCTGCTGGCCGGCACA") and (Rroi[-24:] == "CAGAGCAAGAGAGTGGACTTTTGC") and (len(Rroi[24:-24]) == 462) and ("N" not in Rroi)): # Only include those that have the correct forward primer sequence, correct reverse primer sequence, and the correct number of nucleotides between the primers
        Rroi = Rroi[24:-24] # Trim forward and reverse primers
        variants.append(Rroi)
    #if record_count == 1000: break
  return Counter(variants)

def Output(A107_low_rep1_dict, A107_med_rep1_dict, A107_hi_rep1_dict, A214_low_rep1_dict, A214_med_rep1_dict, A214_hi_rep1_dict, A218_low_rep1_dict, A218_med_rep1_dict,A218_hi_rep1_dict, A107_low_rep2_dict, A107_med_rep2_dict, A107_hi_rep2_dict, A214_low_rep2_dict, A214_med_rep2_dict, A214_hi_rep2_dict, A218_low_rep2_dict, A218_med_rep2_dict,A218_hi_rep2_dict, outfile):
  print ("Compiling results into %s" % outfile)
  outfile = open(outfile,'w')
  muts = list(set(list(A107_low_rep1_dict.keys())+
                  list(A107_med_rep1_dict.keys())+
                  list(A107_hi_rep1_dict.keys())+
                  list(A214_low_rep1_dict.keys())+
                  list(A214_med_rep1_dict.keys())+
                  list(A214_hi_rep1_dict.keys())+
                  list(A218_low_rep1_dict.keys())+
                  list(A218_med_rep1_dict.keys())+
                  list(A218_hi_rep1_dict.keys())+
                  list(A107_low_rep2_dict.keys())+
                  list(A107_med_rep2_dict.keys())+
                  list(A107_hi_rep2_dict.keys())+
                  list(A214_low_rep2_dict.keys())+
                  list(A214_med_rep2_dict.keys())+
                  list(A214_hi_rep2_dict.keys())+
                  list(A218_low_rep2_dict.keys())+
                  list(A218_med_rep2_dict.keys())+
                  list(A218_hi_rep2_dict.keys())))
  outfile.write("\t".join(['mut','A107_low_rep1_count','A107_med_rep1_count','A107_hi_rep1_count', 'A214_low_rep1_count', 'A214_med_rep1_count', 'A214_hi_rep1_count', 'A218_low_rep1_count', 'A218_med_rep1_count', 'A218_hi_rep1_count', 'A107_low_rep2_count', 'A107_med_rep2_count', 'A107_hi_rep2_count', 'A214_low_rep2_count', 'A214_med_rep2_count', 'A214_hi_rep2_count', 'A218_low_rep2_count', 'A218_med_rep2_count', 'A218_hi_rep2_count'])+"\n")
  for mut in muts:
    A107_low_rep1_count   = A107_low_rep1_dict[mut]
    A107_med_rep1_count = A107_med_rep1_dict[mut]
    A107_hi_rep1_count = A107_hi_rep1_dict[mut]
    A214_low_rep1_count = A214_low_rep1_dict[mut]
    A214_med_rep1_count = A214_med_rep1_dict[mut]
    A214_hi_rep1_count = A214_hi_rep1_dict[mut]
    A218_low_rep1_count = A218_low_rep1_dict[mut]
    A218_med_rep1_count = A218_med_rep1_dict[mut]
    A218_hi_rep1_count = A218_hi_rep1_dict[mut]
    A107_low_rep2_count = A107_low_rep2_dict[mut]
    A107_med_rep2_count = A107_med_rep2_dict[mut]
    A107_hi_rep2_count = A107_hi_rep2_dict[mut]
    A214_low_rep2_count = A214_low_rep2_dict[mut]
    A214_med_rep2_count = A214_med_rep2_dict[mut]
    A214_hi_rep2_count = A214_hi_rep2_dict[mut]
    A218_low_rep2_count = A218_low_rep2_dict[mut]
    A218_med_rep2_count = A218_med_rep2_dict[mut]
    A218_hi_rep2_count = A218_hi_rep2_dict[mut]
    outfile.write("\t".join(map(str,[mut,A107_low_rep1_count,A107_med_rep1_count,A107_hi_rep1_count,A214_low_rep1_count,A214_med_rep1_count,A214_hi_rep1_count,A218_low_rep1_count,A218_med_rep1_count,A218_hi_rep1_count,A107_low_rep2_count,A107_med_rep2_count,A107_hi_rep2_count,A214_low_rep2_count,A214_med_rep2_count,A214_hi_rep2_count,A218_low_rep2_count,A218_med_rep2_count,A218_hi_rep2_count]))+"\n")
  outfile.close()

def main():
  outfile = 'result/S2HR1_bind_count_trimmed.tsv'
  A107_low_rep1_dict  = ProcessMultilib('fastq_merged/107_low_rep1_ATCACG_L001_merged_001.assembled.fastq')
  A107_med_rep1_dict  = ProcessMultilib('fastq_merged/107_med_rep1_CGATGT_L001_merged_001.assembled.fastq')
  A107_hi_rep1_dict  = ProcessMultilib('fastq_merged/107_hi_rep1_TTAGGC_L001_merged_001.assembled.fastq')
  A214_low_rep1_dict  = ProcessMultilib('fastq_merged/214_low_rep1_TGACCA_L001_merged_001.assembled.fastq') 
  A214_med_rep1_dict  = ProcessMultilib('fastq_merged/214_med_rep1_ACAGTG_L001_merged_001.assembled.fastq')
  A214_hi_rep1_dict  = ProcessMultilib('fastq_merged/214_hi_rep1_GCCAAT_L001_merged_001.assembled.fastq')
  A218_low_rep1_dict  = ProcessMultilib('fastq_merged/218_low_rep1_CAGATC_L001_merged_001.assembled.fastq')
  A218_med_rep1_dict  = ProcessMultilib('fastq_merged/218_med_rep1_ACTTGA_L001_merged_001.assembled.fastq')
  A218_hi_rep1_dict  = ProcessMultilib('fastq_merged/218_hi_rep1_GATCAG_L001_merged_001.assembled.fastq')
  A107_low_rep2_dict  = ProcessMultilib('fastq_merged/107_low_rep2_TAGCTT_L001_merged_001.assembled.fastq') 
  A107_med_rep2_dict  = ProcessMultilib('fastq_merged/107_med_rep2_GGCTAC_L001_merged_001.assembled.fastq')
  A107_hi_rep2_dict  = ProcessMultilib('fastq_merged/107_hi_rep2_CTTGTA_L001_merged_001.assembled.fastq')
  A214_low_rep2_dict  = ProcessMultilib('fastq_merged/214_low_rep2_AGTCAA_L001_merged_001.assembled.fastq')
  A214_med_rep2_dict  = ProcessMultilib('fastq_merged/214_med_rep2_AGTTCC_L001_merged_001.assembled.fastq')
  A214_hi_rep2_dict  = ProcessMultilib('fastq_merged/214_hi_rep2_ATGTCA_L001_merged_001.assembled.fastq')
  A218_low_rep2_dict  = ProcessMultilib('fastq_merged/218_low_rep2_CCGTCC_L001_merged_001.assembled.fastq') 
  A218_med_rep2_dict  = ProcessMultilib('fastq_merged/218_med_rep2_GTAGAG_L001_merged_001.assembled.fastq') 
  A218_hi_rep2_dict  = ProcessMultilib('fastq_merged/218_hi_rep2_GTCCGC_L001_merged_001.assembled.fastq') 
  Output(A107_low_rep1_dict, A107_med_rep1_dict, A107_hi_rep1_dict, A214_low_rep1_dict, A214_med_rep1_dict, A214_hi_rep1_dict, A218_low_rep1_dict, A218_med_rep1_dict,A218_hi_rep1_dict, A107_low_rep2_dict, A107_med_rep2_dict, A107_hi_rep2_dict, A214_low_rep2_dict, A214_med_rep2_dict, A214_hi_rep2_dict, A218_low_rep2_dict, A218_med_rep2_dict,A218_hi_rep2_dict, outfile)

if __name__ == "__main__":
  main()
