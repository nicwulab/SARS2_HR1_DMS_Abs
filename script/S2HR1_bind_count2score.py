#!/usr/bin/python
import sys
import pandas as pd
import numpy as np

def count_to_freq(df, colname):
    df[colname+'_freq'] = (df[colname]+1)/(df[colname].sum()+len(df))
    return (df)

def A107_score_calculate(df, rep):
    print (rep)
    A107_weight = df['A107_low_'+rep+'_count_freq']*0.33 + df['A107_med_'+rep+'_count_freq']*0.67 + \
                 df['A107_hi_'+rep+'_count_freq']*1
    A107_norm_factor = df['A107_low_'+rep+'_count_freq'] + df['A107_med_'+rep+'_count_freq'] + \
                 df['A107_hi_'+rep+'_count_freq']
    A107_avg_freq = A107_norm_factor/4
    df['A107_avg_freq_'+rep] = A107_avg_freq
    df['A107_weight_'+rep] = A107_weight/A107_norm_factor
    w_summary = df.groupby('mut_class')['A107_weight_'+rep].mean()
    w_summary = w_summary.reset_index()
    print (w_summary)
    w_silent   = (float(w_summary.loc[w_summary['mut_class']=='silent']['A107_weight_'+rep]))
    w_nonsense = (float(w_summary.loc[w_summary['mut_class']=='nonsense']['A107_weight_'+rep]))
    df['A107_score_'+rep] = (df['A107_weight_'+rep]-w_nonsense)/(w_silent-w_nonsense)
    print ('w_nonsense', w_nonsense)
    print ('w_silent', w_silent)
    return (df)

def A214_score_calculate(df, rep):
    print (rep)
    A214_weight = df['A214_low_'+rep+'_count_freq']*0.33 + df['A214_med_'+rep+'_count_freq']*0.67 + \
                 df['A214_hi_'+rep+'_count_freq']*1
    A214_norm_factor = df['A214_low_'+rep+'_count_freq'] + df['A214_med_'+rep+'_count_freq'] + \
                 df['A214_hi_'+rep+'_count_freq']
    A214_avg_freq = A214_norm_factor/4
    df['A214_avg_freq_'+rep] = A214_avg_freq
    df['A214_weight_'+rep] = A214_weight/A214_norm_factor
    w_summary = df.groupby('mut_class')['A214_weight_'+rep].mean()
    w_summary = w_summary.reset_index()
    print (w_summary)
    w_silent   = (float(w_summary.loc[w_summary['mut_class']=='silent']['A214_weight_'+rep]))
    w_nonsense = (float(w_summary.loc[w_summary['mut_class']=='nonsense']['A214_weight_'+rep]))
    df['A214_score_'+rep] = (df['A214_weight_'+rep]-w_nonsense)/(w_silent-w_nonsense)
    print ('w_nonsense', w_nonsense)
    print ('w_silent', w_silent)
    return (df)

def A218_score_calculate(df, rep):
    print (rep)
    A218_weight = df['A218_low_'+rep+'_count_freq']*0.33 + df['A218_med_'+rep+'_count_freq']*0.67 + \
                 df['A218_hi_'+rep+'_count_freq']*1
    A218_norm_factor = df['A218_low_'+rep+'_count_freq'] + df['A218_med_'+rep+'_count_freq'] + \
                 df['A218_hi_'+rep+'_count_freq']
    A218_avg_freq = A218_norm_factor/4
    df['A218_avg_freq_'+rep] = A218_avg_freq
    df['A218_weight_'+rep] = A218_weight/A218_norm_factor
    w_summary = df.groupby('mut_class')['A218_weight_'+rep].mean()
    w_summary = w_summary.reset_index()
    print (w_summary)
    w_silent   = (float(w_summary.loc[w_summary['mut_class']=='silent']['A218_weight_'+rep]))
    w_nonsense = (float(w_summary.loc[w_summary['mut_class']=='nonsense']['A218_weight_'+rep]))
    df['A218_score_'+rep] = (df['A218_weight_'+rep]-w_nonsense)/(w_silent-w_nonsense)
    print ('w_nonsense', w_nonsense)
    print ('w_silent', w_silent)
    return (df)

def read_count_data(count_file):
    print ('reading: %s' % count_file)
    df = pd.read_csv(count_file, sep='\t')
    colnames = [colname for colname in df]
    for colname in colnames:
        if 'mut' not in colname and 'frag' not in colname:
            df = count_to_freq(df, colname)
    df = A107_score_calculate(df, 'rep1')
    df = A107_score_calculate(df, 'rep2')
    df = A214_score_calculate(df, 'rep1')
    df = A214_score_calculate(df, 'rep2')
    df = A218_score_calculate(df, 'rep1')
    df = A218_score_calculate(df, 'rep2')
    df['A107_score'] = (df['A107_score_rep1'] + df['A107_score_rep2'])/2
    df['A214_score'] = (df['A214_score_rep1'] + df['A214_score_rep2'])/2
    df['A218_score'] = (df['A218_score_rep1'] + df['A218_score_rep2'])/2
    df['avg_freq'] = (df['A107_avg_freq_rep1'] + df['A107_avg_freq_rep2'] + \
                      df['A214_avg_freq_rep1'] + df['A214_avg_freq_rep2'] + \
                      df['A218_avg_freq_rep1'] + df['A218_avg_freq_rep2'])/6
    return (df)

def main():
    outfile_1 = "result/S2HR1_bind_scores.tsv"
    df = read_count_data("result/S2HR1_bind_count_trimmed_aa.tsv")
    print ('writing: %s' % outfile_1)
    df.to_csv(outfile_1, sep="\t", index=False)

if __name__ == "__main__":
    main()
