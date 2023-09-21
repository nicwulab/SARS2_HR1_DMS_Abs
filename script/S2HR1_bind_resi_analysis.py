#!/usr/bin/python3
import numpy as np
import pandas as pd

def rms(values):
  return np.sqrt(sum(values**2)/len(values))

def main():
  freq_cutoff = 0.00002
  outfile = 'result/S2HR1_adj_score_by_resi.tsv'
  df = pd.read_csv('result/S2HR1_scores_common.tsv',sep="\t")
  df = df[df['avg_freq'] >= freq_cutoff]
  df = df[df['mut_class'] == 'missense']
  df['resi'] = df['mut'].str[0:-1]
  df_resi_mean  = df.groupby('resi')['residual'].mean().reset_index(name='mean_adjust_bind')
  df_resi_mean['pos'] = df_resi_mean['resi'].str[1::].astype(int)
  df_resi_mean = df_resi_mean.sort_values(by=['pos'])
  print ('writing: %s' % outfile)
  df_resi_mean.to_csv(outfile, sep="\t", index=False)

if __name__ == "__main__":
  main()
