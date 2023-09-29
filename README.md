## Characterizing the epitope of a public antibody clonotype to SARS-CoV-2 spike HR1 using deep mutational scanning
This README describes the analysis in:   
[Evidence of antigenic drift in the fusion machinery core of SARS-CoV-2 spike](https://www.biorxiv.org/content/10.1101/2023.09.27.559757v1)

### Input files:
* Raw reads from [BioProject PRJNA888135](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA888135)
* [./data/barcodes.tsv](./data/barcodes.tsv): Internal barcodes for NGS error correction
* [./Fasta/HR1_ref_seq.fa](./Fasta/HR1_ref_seq.fa): Reference (i.e. wild type) amino acid sequences
* [./data/S2HR1_exp_fus_scores.tsv](./data/S2HR1_exp_fus_scores.tsv): Expression and fusion scores from [a previous study](https://www.nature.com/articles/s41467-023-37786-1)

### Computing binding scores from raw Illumina MiSeq reads
1. Use [PEAR](https://github.com/tseemann/PEAR) to merge forward and reverse reads.
2. Run ``python script/S2HR1_bind_fastq2count.py`` to generate sequence count file.
3. Run ``python script/S2HR1_bind_count_nuc2aa.py`` to generate [amino acid count file](./result/S2HR1_bind_count_trimmed_aa.tsv). 
4. Run ``python script/S2HR1_bind_count2score.py`` to generate [binding score file](./result/S2HR1_bind_scores.tsv).
5. Run ``Rscript script/common_muts.R`` to add expression and fusion scores into the [final score file](./result/S2HR1_scores_common.tsv).
6. Run ``python script/S2HR1_bind_resi_analysis.py`` to generate [mutational effect file](./result/S2HR1_adj_score_by_resi.tsv)
7. Run ``python script/Bfactor_to_score.py`` to generate [PDB with mutational effect as B factor](PDB/6vxx_score.pdb)

### Plot correlations between replicates and between antibodies
1. Run ``Rscript script/plot_correlations.R`` to compare binding scores across antibodies and replicates.
2. Run ``Rscript script/plot_QC.R`` to compare binding scores among nonsense, missense, and silent mutations.
3. Run ``Rscript script/plot_bind_vs_exp`` to compare binding and expression scores.

### Analyzing mutational effects 
1. Run ``Rscript script/plot_adj_bind_score.R`` to plot mutational effects along residue positions.
2. Run ``pymol script/plot_Bfactor_as_score.pml`` to generate structural figure colored according to mutational effects.

### Expeirmental validation results
* See [./flow_cytometry](./flow_cytometry)
