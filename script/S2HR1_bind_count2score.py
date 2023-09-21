#!/usr/bin/python
import sys
import pandas as pd
import numpy as np

def count_to_freq(df, colname):
    df[colname+'_freq'] = (df[colname]+1)/(df[colname].sum()+len(df))
    return (df)

def cell_counting(s):
  if s=='low': return (100)
  if s=='med': return (10)
  if s=='hi':  return (1)

def MFI(s):
  if s=='low':  return (1)
  if s=='med':  return (50)
  if s=='hi':  return (500)

def score_calculate(df, Ab, rep, freq_cutoff):
    print (rep)
    weight     = df[Ab+'_low_'+rep+'_count_freq']*cell_counting('low')*MFI('low') + \
                 df[Ab+'_med_'+rep+'_count_freq']*cell_counting('med')*MFI('med') + \
                 df[Ab+'_hi_'+rep+'_count_freq']*cell_counting('hi')*MFI('hi')
    cell_count = df[Ab+'_low_'+rep+'_count_freq']*cell_counting('low') + \
                 df[Ab+'_med_'+rep+'_count_freq']*cell_counting('med') + \
                 df[Ab+'_hi_'+rep+'_count_freq']*cell_counting('hi')
    df[Ab+'_avg_freq_'+rep] = cell_count/(cell_counting('low')+cell_counting('med')+cell_counting('hi'))
    df[Ab+'_weight_'+rep] = np.log10(weight/cell_count)
    df_high_freq = df[df[Ab+'_avg_freq_'+rep] >= freq_cutoff]
    w_summary = df_high_freq.groupby('mut_class')[Ab+'_weight_'+rep].mean()
    w_summary = w_summary.reset_index()
    print (w_summary)
    w_silent   = (float(w_summary.loc[w_summary['mut_class']=='silent'][Ab+'_weight_'+rep]))
    w_nonsense = (float(w_summary.loc[w_summary['mut_class']=='nonsense'][Ab+'_weight_'+rep]))
    df[Ab+'_score_'+rep] = (df[Ab+'_weight_'+rep]-w_nonsense)/(w_silent-w_nonsense)
    print ('w_nonsense', w_nonsense)
    print ('w_silent', w_silent)
    return (df)

def read_count_data(count_file, freq_cutoff):
    print ('reading: %s' % count_file)
    df = pd.read_csv(count_file, sep='\t')
    colnames = [colname for colname in df]
    for colname in colnames:
        if 'mut' not in colname and 'frag' not in colname:
            df = count_to_freq(df, colname)
    df = score_calculate(df, 'A107', 'rep1', freq_cutoff)
    df = score_calculate(df, 'A107', 'rep2', freq_cutoff)
    df = score_calculate(df, 'A214', 'rep1', freq_cutoff)
    df = score_calculate(df, 'A214', 'rep2', freq_cutoff)
    df = score_calculate(df, 'A218', 'rep1', freq_cutoff)
    df = score_calculate(df, 'A218', 'rep2', freq_cutoff)
    df['A107_score'] = (df['A107_score_rep1'] + df['A107_score_rep2'])/2
    df['A214_score'] = (df['A214_score_rep1'] + df['A214_score_rep2'])/2
    df['A218_score'] = (df['A218_score_rep1'] + df['A218_score_rep2'])/2
    df['avg_freq'] = (df['A107_avg_freq_rep1'] + df['A107_avg_freq_rep2'] + \
                      df['A214_avg_freq_rep1'] + df['A214_avg_freq_rep2'] + \
                      df['A218_avg_freq_rep1'] + df['A218_avg_freq_rep2'])/6
    return (df)

def main():
    outfile_1 = "result/S2HR1_bind_scores.tsv"
    outfile_2 = "result/S2HR1_mean_bind_scores.tsv"
    freq_cutoff = 0.00002
    df = read_count_data("result/S2HR1_bind_count_trimmed_aa.tsv", freq_cutoff)
    print ('writing: %s' % outfile_1)
    df.to_csv(outfile_1, sep="\t", index=False)
    df['resi'] = df['mut'].str[0:-1]
    all_resi   = df[df['mut_class'] != 'WT'].groupby('resi').size().reset_index()
    df_by_resi = df
    df_by_resi = df_by_resi[df_by_resi['avg_freq'] >= freq_cutoff]
    print ('# of missense:', len(df_by_resi[df_by_resi['mut_class'] == 'missense']))
    print ('# of silent:', len(df_by_resi[df_by_resi['mut_class'] == 'silent']))
    print ('# of nonsense:', len(df_by_resi[df_by_resi['mut_class'] == 'nonsense']))
    df_by_resi = df_by_resi[df_by_resi['mut_class'] != 'WT']
    df_by_resi = df_by_resi[df_by_resi['mut_class'] != 'silent']
    df_by_resi = df_by_resi[df_by_resi['mut_class'] != 'nonsense']
    df_by_resi_mean_A107  = df_by_resi.groupby('resi')['A107_score'].mean().reset_index(name='mean_A107_score')
    df_by_resi_mean_A214  = df_by_resi.groupby('resi')['A214_score'].mean().reset_index(name='mean_A214_score')
    df_by_resi_mean_A218  = df_by_resi.groupby('resi')['A218_score'].mean().reset_index(name='mean_A218_score')
    df_by_resi_mean = pd.merge(df_by_resi_mean_A107, df_by_resi_mean_A214, on='resi', how='outer')
    df_by_resi_mean = pd.merge(df_by_resi_mean, df_by_resi_mean_A218, on='resi', how='outer')
    df_by_resi_count = df_by_resi.groupby('resi').size().reset_index(name='count')
    df_by_resi = pd.merge(df_by_resi_mean, df_by_resi_count, on='resi', how='outer')
    df_by_resi = pd.merge(all_resi, df_by_resi, on='resi', how='outer')
    df_by_resi = df_by_resi.sort_values(by='resi', key=lambda x:x.str[1::].astype(int))
    df_by_resi['pos'] = df_by_resi['resi'].str[1::].astype(int)
    df_by_resi = df_by_resi[['resi','pos','count','mean_A107_score','mean_A214_score','mean_A218_score']]
    print ('writing: %s' % outfile_2)
    df_by_resi.to_csv(outfile_2, sep="\t", index=False)

if __name__ == "__main__":
    main()
