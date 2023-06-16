#!/usr/bin/python
import os
import sys
import numpy as np
import pandas as pd
from functools import partial
from collections import defaultdict

def translation(seq):
  dnamap = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"_", "TAG":"_",
    "TGT":"C", "TGC":"C", "TGA":"_", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",
    "NNK":"X"}
  pep = []
  i = 0
  while i < len(seq):
    codon = seq[i:i+3]
    aa = dnamap[codon]
    pep.append(aa)
    i = i + 3
  pep = ''.join(pep)
  return pep

def reading_barcode(barfile):
  print ('reading: %s' % barfile)
  infile = open(barfile, 'r')
  bc_dict = {}
  for line in infile.readlines():
    if 'barcode' in line: continue
    bc, ID = line.rstrip().rsplit("\t")
    if ID == 'WT': 
      bc_dict[bc] = 'WT'
    else:
      cas, cas_pos = list(map(int,ID.replace('Cassette','').rsplit('_')))
      pos = (cas-1)*8+cas_pos
      bc_dict[bc] = pos
  return (bc_dict)

def nuc_2_aa(countfile, bc_dict, refseq):
  print ('reading: %s' % countfile)
  infile    = open(countfile, 'r')
  count_line = 0
  data_dict  = {}
  for line in infile.readlines():
    count_line += 1
    if count_line == 1:
      header = line.rstrip()
      continue
    line = line.rstrip().rsplit("\t")
    DNA  = line[0][3:-3]
    data = np.array(list(map(int,line[1::])))
    if len(DNA) != 456: continue
    barcode = DNA[2::3] 
    if barcode not in bc_dict.keys(): continue
    mut_pos = bc_dict[barcode]
    pep = translation(DNA)
    mut = 'WT' if mut_pos == 'WT' else refseq[int(mut_pos)-1]+str(mut_pos)+pep[int(mut_pos)-1]
    if mut not in data_dict.keys():
      data_dict[mut] = data
    else:
      data_dict[mut] += data
  infile.close()
  return header, data_dict

def mut_classification(mut):
  if mut=='WT':
    return 'WT'
  elif mut[-1] == '_':
    return 'nonsense'
  elif mut[0]==mut[-1]: 
    return 'silent'
  else:
    return 'missense'

def writing_aa_count(data_dict, header, outfile):
  print ('writing: %s' % outfile)
  outfile   = open(outfile, 'w') 
  header = ['mut','mut_class']+header.rsplit("\t")[1::]
  outfile.write("\t".join(header)+"\n")
  muts = list(data_dict.keys())
  muts.remove('WT')
  for mut in sorted(muts, key=lambda x:int(x[1:-1]))+['WT']:
    adj_mut   = mut if 'WT' in mut else mut[0]+str(int(mut[1:-1])+882)+mut[-1]
    mut_class = mut_classification(adj_mut)
    outfile.write("\t".join(map(str, [adj_mut,mut_class]))+"\t"+
                  "\t".join(list(map(str,data_dict[mut])))+"\n")
  outfile.close()

def main():
  outfile   = 'result/S2HR1_bind_count_trimmed_aa.tsv'
  countfile = 'result/S2HR1_bind_count_trimmed.tsv'
  barfile   = 'data/barcodes.tsv'
  refseq    = [line for line in open('Fasta/HR1_ref_seq.fa','r').readlines()][1].rstrip()
  bc_dict   = reading_barcode(barfile)
  header, data_dict = nuc_2_aa(countfile, bc_dict, refseq)
  writing_aa_count(data_dict, header, outfile)

if __name__ == "__main__":
  main()
