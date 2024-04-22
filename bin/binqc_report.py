#!/usr/bin/env python

import os
import sys
import numpy as np
import pandas as pd


dbname = sys.argv[1]
checkm_summary_file = sys.argv[2]
gtdbtk_summary_file = sys.argv[3]
genome_summary_file = sys.argv[4]
rrna_summary_file = sys.argv[5]
gunc_summary_file = sys.argv[6]

# Output file name
report_csv = f'{dbname}.report.csv'

def checkm_score(row):
    # From GTDB:
    # Completeness - (5 * Contamination)
    # Should be less than or equal to 50
    return row['Completeness'] - (5 * row['Contamination'])

## GTDB
def classify_method_stats(row):
    row['bin_name'] = row['user_genome']
    if row['classification_method']=="ANI":
        row['ani']=row['fastani_ani']
        row['aln_frac']=row['fastani_af']
        row['closest_ref']=row['fastani_reference']
    elif row['classification_method']=="ANI/Placement":
        row['ani']=row['closest_placement_ani']
        row['aln_frac']=row['closest_placement_af']
        row['closest_ref']=row['closest_placement_reference']
    else:
        row['ani']=row['fastani_ani']
        row['aln_frac']=row['fastani_af']
        row['closest_ref']=row['fastani_reference']
    return row

## CheckM
def parse_checkm_output(checkm_summary):
    df = pd.read_table(checkm_summary, comment='[').rename(columns={"Bin Id": "bin_name", "Strain heterogeneity":"Strain_heterogeneity"})
    df['checkm_score'] = df.apply(checkm_score, axis = 1)
    keep_cols = ['bin_name','Completeness', 'Contamination', 'Strain_heterogeneity', 'checkm_score']
    return df[keep_cols].sort_values('checkm_score', ascending = True).reset_index(drop = True)

# GTDB
def parse_gtdb_outputs(gtdb_output):
    df = pd.read_table(gtdb_output)
    df = df.apply(classify_method_stats, axis = 1)
    keep_col = ['bin_name','classification','classification_method','closest_ref', 'ani', 'aln_frac','msa_percent','red_value', 'note', 'warnings']
    return df[keep_col].set_index('bin_name')

# GUNC
def parse_gunc_outputs(gunc_output):
    df = pd.read_table(gunc_output).rename(columns={"genome":"bin_name", "pass.GUNC":"pass_GUNC"} )
    keep_col = ['bin_name','pass_GUNC']
    return df[keep_col].set_index('bin_name')

#rRna
def parse_rrna_stats(genes_file):
    return pd.read_table(genes_file, names=['bin_name','num_rRNAs'],
                        index_col='bin_name')

# SeqKit
def parse_seqkit_stats(stats_file):
    df = pd.read_table(stats_file, names=['file', 'format', 'type', 'num_seqs', 'sum_len', 'min_len', 'avg_len', 'max_len'])
    df['bin_name'] = df['file'].apply(lambda x: os.path.basename(os.path.splitext(x)[0]))
    keep_cols=['bin_name','num_seqs','sum_len','min_len','avg_len','max_len']
    return df[keep_cols].set_index('bin_name')

# Quast
def parse_quast_stats(quast_file):
    df = pd.read_table(quast_file).rename(columns={"Assembly": "bin_name", "# predicted rRNA genes":"#_predicted_rRNA_genes", "GC (%)": "GC(%)" })
    df['bin_name'] = df['bin_name'].apply(lambda x: x.replace("_", "-"))
    keep_cols=['bin_name','#_predicted_rRNA_genes', 'GC(%)']
    return df[keep_cols].set_index('bin_name')

# Quast-Longest Contig
def parse_quast_stats2(quast_file):
    df = pd.read_table(quast_file).rename(columns={"Assembly": "bin_name", "Largest contig":"longest_contig"})
    df['bin_name'] = df['bin_name'].apply(lambda x: x.replace("_", "-"))
    keep_cols=['bin_name','longest_contig']
    return df[keep_cols].set_index('bin_name')

if __name__ == '__main__':
    print('Parsing CheckM output ...')
    checkm_df = parse_checkm_output(checkm_summary_file).set_index('bin_name')

    print('Parsing GTDBtk output ...')
    gtdb_df = parse_gtdb_outputs(gtdbtk_summary_file)
    
    print('Parsing Seqkit stats ...')
    bin_stats_df = parse_seqkit_stats(genome_summary_file)

    print('Parsing Barrnap output ...')
    rrna_df = parse_rrna_stats(rrna_summary_file)
    
    print('Parsing GUNC output ...')
    gunc_df = parse_gunc_outputs(gunc_summary_file)

    print(rrna_df)

    print(checkm_df)

    print(gtdb_df)

    print(bin_stats_df)
    ## Aggregate
    print('Aggregating data ...')
    bins_df = pd.concat([bin_stats_df,
                         gtdb_df,
                         checkm_df,
                         rrna_df,
                         gunc_df],
                         axis=1, sort=False).sort_values('checkm_score', ascending = True).reset_index().rename(columns={'index':'bin_name'})
    print('Writing output ...')
    bins_df.to_csv(report_csv, index = False)
    print('Done.')
