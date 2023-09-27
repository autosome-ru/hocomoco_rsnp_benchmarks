#!/usr/bin/python3.8
import os
import sys
import subprocess
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, kendalltau, fisher_exact


BASE_PATH = '/'.join(os.getcwd().split('/')[:-1])
REF_PATH = BASE_PATH + '/assembly/hg38.fa'
PWM_PATH = BASE_PATH + \
    '/motifs/pwm/{name}/{name}@{pwm_name}.pwm'
THR_PATH = BASE_PATH + \
    '/motifs/thr/{name}@{pwm_name}.thr'
FASTA_PATH = BASE_PATH + \
    '/fasta/{name}@{pwm_name}@a.fasta'
PV_PATH = BASE_PATH + \
    '/fasta/{name}@{pwm_name}@a.processed.fasta'
ADASTRA_PATH = BASE_PATH + '/adastra/TF'


TP_THR = 0.05  # True positive threshold for min FDR-adjusted pvalue
TN_THR = 0.5  # True negative threshold for min FDR-adjusted pvalue
LH_THR = -np.log10(0.0005)  # Likelihood threshold for motif entry


def NOT_FOUND(name, pwm_name):
    return {
        'name': name,
        'pwm': pwm_name,
        'total': 0,
        'positive': 0,
        'negative': 0,
        'concordant': 0,
        'discordant': 0,
        'odds_ratio': None,
        'fisher_p': None,
        'r': None,
        'tau': None,
    }


def main(name, pwm_name):
    df = parse_adastra(name)
    if df.empty:
        return NOT_FOUND(name, pwm_name)
    df = filter_by_annotation(df)
    df = df.apply(lambda x: add_ground_truth(x), axis=1)
    if len(df.index) == 0:
        return NOT_FOUND(name, pwm_name)
    df.dropna(subset=['true'], inplace=True)
    if len(df.index) == 0:
        return NOT_FOUND(name, pwm_name)
    pwm_path = PWM_PATH.format(name=name, pwm_name=pwm_name)
    thr_path = THR_PATH.format(name=name, pwm_name=pwm_name)
    motif_length = get_motif_length(pwm_path)
    fasta_path = FASTA_PATH.format(name=name, pwm_name=pwm_name)
    pval_path = PV_PATH.format(name=name, pwm_name=pwm_name)
    make_fasta_file(df, fasta_path, motif_length)
    make_sarus_file(fasta_path, pval_path, pwm_path, thr_path)
    snp_dict = get_snp_dict(pval_path)
    df = df.apply(lambda x: add_scores(x, snp_dict), axis=1)
    total = len(df.index)
    positive = df[df['true'] == 1]
    negative = df[df['true'] == 0]
    contingency_table = get_contingency_table(positive, negative)
    df.dropna(subset=['logpbsp'], inplace=True)
    logfc = df['logfc'].to_list()
    logpbsp = df['logpbsp'].to_list()
    fisher_test = fisher_exact(contingency_table)
    r = None
    tau = None
    if len(logfc) > 2 and len(logpbsp) > 2:
        r = pearsonr(logfc, logpbsp)[0]
        tau = kendalltau(logfc, logpbsp)[0]
    return {
        'name': name,
        'pwm': pwm_name,
        'total': total,
        'positive': len(positive.index),
        'negative': len(negative.index),
        'odds_ratio': fisher_test[0],
        'fisher_p': fisher_test[1],
        'r': r,
        'tau': tau,
    }


def parse_adastra(name):
    df = pd.DataFrame()
    try:
        df = pd.read_csv(ADASTRA_PATH + f'/{name}_HUMAN.tsv', sep='\t')
    except FileNotFoundError:
        return df
    return df


def filter_by_annotation(df):
    if df.empty:
        return df
    return df[(df['motif_log_pref'] > LH_THR) |
              (df['motif_log_palt'] > LH_THR)]


def add_ground_truth(row):
    if min(row['fdrp_bh_ref'], row['fdrp_bh_alt']) < TP_THR:
        row['true'] = True
    elif min(row['fdrp_bh_ref'], row['fdrp_bh_alt']) > TN_THR:
        row['true'] = False
    else:
        row['true'] = None
    return row


def get_motif_length(pwm_path):
    with open(pwm_path, 'r') as f:
        motif_length = len(f.readlines()) - 1
        if motif_length == 0 or motif_length is None:
            print('Invalid PWM')
            sys.exit(0)
        return motif_length


def make_fasta_file(df, fasta, motif_length):
    if not os.path.isfile(fasta):
        if motif_length is None:
            print('No such PWM')
            sys.exit(0)
        with open(fasta, 'w') as outfile:
            df = df.apply(
                lambda row: write_sequence(row, motif_length, outfile), axis=1)


def write_sequence(row, motif_length, outfile):
    chr = row['#chr']
    pos = row['pos']
    ref = row['ref']
    alt = row['alt']
    start = pos - motif_length
    end = pos + motif_length - 1
    id = f'{row["ID"]}_{row["fdrp_bh_alt"]:.5f}'
    coords = subprocess.Popen(['echo', f'{chr}\t{start}\t{end}'],
                              stdout=subprocess.PIPE)
    output = subprocess.check_output(['../external_programs/bedtools.static',
                                      'getfasta', '-fi', REF_PATH, '-bed',
                                      'stdin'],
                                     stdin=coords.stdout)
    sequence = str(output, 'utf-8').split('\n')[1]
    if ref.upper() != sequence[motif_length - 1].upper():
        print(f'{id}: Failed')
    else:
        left = sequence[:motif_length-1].upper()
        right = sequence[motif_length:].upper()
        outfile.write(f'>{id}_ref\n{left}{ref}{right}\n')
        outfile.write(f'>{id}_alt\n{left}{alt}{right}\n')


def make_sarus_file(fasta, pvalues, pwm, thr):
    if not os.path.isfile(pvalues):
        output = get_pvalues(fasta, pwm, thr)
        with open(pvalues, 'wb') as out:
            out.write(output)


def get_pvalues(fasta_path, pwm_path, thr_path):
    return subprocess.check_output(['java', '-cp',
                                    '../external_programs/sarus-2.0.2.jar',
                                    'ru.autosome.SARUS',
                                    fasta_path, pwm_path,
                                    'besthit',
                                    '--pvalues-file', thr_path,
                                    '--threshold-mode', 'score',
                                    '--output-scoring-mode',
                                    'logpvalue'])


def get_snp_dict(pvalues_path):
    snp_dict = {}
    if os.path.isfile(pvalues_path) is not None:
        with open(pvalues_path, 'r') as sarus:
            allele = None
            current_snp_id = None
            for line in sarus:
                if line[0] == '>':
                    allele = line[-4:-1]
                    current_snp_id = line[1:-5]
                    assert allele in ('ref', 'alt')
                    if allele == 'ref':
                        snp_dict[current_snp_id] = \
                            {'ref': {'p': None, 'pos': -1, 'orient': None},
                             'alt': {'p': None, 'pos': -1, 'orient': None}}
                else:
                    assert allele in ('ref', 'alt')
                    line = line.strip('\n').split()
                    snp_dict[current_snp_id][allele] = {
                        'p': float(line[0]),
                        'pos': int(line[1]),
                        'orient': line[2]
                    }
    return snp_dict


def add_scores(row, snp_dict):
    ID = f'{row["ID"]}_{row["fdrp_bh_alt"]:.5f}'
    ref_p = snp_dict[ID]['ref']['p']
    alt_p = snp_dict[ID]['alt']['p']
    best_p = max(ref_p, alt_p)
    row['best_p'] = best_p
    row['logfc'] = alt_p - ref_p
    row['score'] = abs(row['logfc'])
    row['motif_conc'] = np.sign((row['fdrp_bh_ref'] - row['fdrp_bh_alt']) *
                                np.sign(row['logfc']))
    if min(row['fdrp_bh_ref'], row['fdrp_bh_alt']) != 0:
        row['logpbsp'] = \
            -np.log10(min(row['fdrp_bh_ref'], row['fdrp_bh_alt'])) * \
            np.sign((row['fdrp_bh_ref'] - row['fdrp_bh_alt']))
    else:
        row['logpbsp'] = None
    return row


def get_contingency_table(positive, negative):
    tp = len(positive[positive['motif_conc'] == 1].index)
    fp = len(positive[positive['motif_conc'] == -1].index)
    tn = len(negative[negative['motif_conc'] == 1].index)
    fn = len(negative[negative['motif_conc'] == -1].index)
    return [[tp, fp], [fn, tn]]


if __name__ == '__main__':
    with open('adastra_motifs.tsv', 'a') as results:
        tf_name = sys.argv[1].split('@')[0]
        pwm_name = '@'.join(sys.argv[1].split('@')[1:])
        result = main(tf_name, pwm_name)
        row = [result['name'],
               result['pwm'],
               str(result['total']),
               str(result['positive']),
               str(result['negative']),
               str(result['odds_ratio']),
               str(result['fisher_p']),
               str(result['r']),
               str(result['tau'])]
        results.write('\t'.join(row) + '\n')
