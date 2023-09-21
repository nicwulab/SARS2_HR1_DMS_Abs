#!/usr/bin/python
import os
import sys
from collections import defaultdict

def read_scorefile(file_file):
  infile = open(file_file, 'r')
  score_dict = {}
  for line in infile.readlines():
    if 'pos' in line: continue
    resi, score, pos = line.rstrip().rsplit("\t")
    score_dict[pos] = float(score)
  infile.close()
  return (score_dict)

def normalizing_score(score_dict):
  max_score  = max(score_dict.values())
  min_score  = min(score_dict.values())
  print ('score range (pre-norm): %f to %f' % (min_score, max_score))
  norm_score = []
  for pos in score_dict.keys():
    score = score_dict[pos]
    score = (score-min_score)/(max_score-min_score)*100
    score_dict[pos] = score
    norm_score.append(score)
  print ('score range (post-norm): %f to %f' % (min(norm_score), max(norm_score)))
  return (score_dict)

def add_score_to_pdb(score_dict, pdb_file):
  assert('.pdb' in pdb_file)
  new_pdb_file = pdb_file.replace('.pdb', '_score.pdb')
  print ("writing: %s" % new_pdb_file)
  infile  = open(pdb_file, 'r')
  outfile = open(new_pdb_file,'w')
  for line in infile.readlines():
    if "ATOM" == line[0:4]:
      pos      = int(line[22:26])
      chain    = line[21:22].rstrip()
      b_factor = line[60:66]
      if chain == 'A' and str(pos) in score_dict.keys():
        score = str(round(score_dict[str(pos)],2))
        score = score+'0'*(2-len(score.rsplit('.')[-1]))
        score = ' '*(6-len(score))+score
        new_line = line[0:60]+score+line[66::]
        outfile.write(new_line)
      else:
        outfile.write(line)
    else:
      outfile.write(line)
  outfile.close()
  
def main():
  pdb_file = "PDB/6vxx.pdb"
  score_file = 'result/S2HR1_adj_score_by_resi.tsv'
  score_dict = read_scorefile(score_file)
  score_dict = normalizing_score(score_dict)
  add_score_to_pdb(score_dict, pdb_file)

if __name__ == "__main__":
  main() 
