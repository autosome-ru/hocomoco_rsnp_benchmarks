#!/usr/bin/python3.8
import os
import sys
import subprocess
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, kendalltau
from sklearn.metrics import roc_auc_score, precision_recall_curve, auc


BASE_PATH = '/'.join(os.getcwd().split('/')[:-1])
REF_PATH = BASE_PATH + '/assembly/hg19.fa'
PWM_PATH = BASE_PATH + \
    '/motifs/pwm/{name}/{name}@{pwm_name}.pwm'
THR_PATH = BASE_PATH + \
    '/greco/thr/{name}@{pwm_name}.thr'
FASTA_PATH = BASE_PATH + \
    '/fasta/{name}@{pwm_name}{prefix}.fasta'
PV_PATH = BASE_PATH + \
    '/fasta/{name}@{pwm_name}{prefix}.processed.fasta'
BATCH1_PATH = BASE_PATH + '/selex/batch1/'
BATCH2_PATH = BASE_PATH + '/selex/batch2/'

TP_THR = 0.01  # True positive threshold
TN_THR = 0.5  # True negative threshold
LH_THR = 0.05  # Likelihood threshold for the 1st batch


def NOT_FOUND(name, pwm_name, batch):
    return {
        'name': name,
        'pwm': pwm_name,
        'batch': batch,
        'total': 0,
        'positive': 0,
        'negative': 0,
        'auroc': None,
        'auprc': None,
        'r': None,
        'tau': None,
        'metric': None
    }


def main(name, pwm_name, batch=1, metric='shift-max'):
    df = parse_snpselex(name, batch)
    if df.empty:
        return NOT_FOUND(name, pwm_name, batch)
    df = df.apply(lambda x: add_ground_truth(x, batch), axis=1)
    df.dropna(subset=['true'], inplace=True)
    if len(df.index) == 0:
        return NOT_FOUND(name, pwm_name)
    pwm_path = PWM_PATH.format(name=name, pwm_name=pwm_name)
    thr_path = THR_PATH.format(name=name, pwm_name=pwm_name)
    motif_length = get_motif_length(pwm_path)
    prefix = '@1' if batch == 1 else ''
    fasta_path = FASTA_PATH.format(name=name, pwm_name=pwm_name, prefix=prefix)
    pval_path = PV_PATH.format(name=name, pwm_name=pwm_name, prefix=prefix)
    make_fasta_file(df, fasta_path, motif_length)
    make_sarus_file(fasta_path, pval_path, pwm_path, thr_path)
    snp_dict = get_snp_dict(pval_path, motif_length)
    df = df.apply(lambda x: add_score_best_p(x, snp_dict, metric), axis=1)
    df.reset_index(inplace=True, drop=True)
    total = len(df.index)
    positive = len(df[df["true"] == 1].index)
    negative = len(df[df["true"] == 0].index)
    df.dropna(subset=['logpbsp'], inplace=True)
    positive_entries = df[df["true"] == 1]
    logfc = positive_entries['logfc'].to_list()
    logpbsp = positive_entries['logpbsp'].to_list()
    true = df['true'].to_list()
    score = df['score'].to_list()
    auroc = None
    auprc = None
    if True in true and False in true:
        auroc = roc_auc_score(true, score)
        precision, recall, thr = precision_recall_curve(true, score)
        auprc = auc(recall, precision)
    r = pearsonr(logfc, logpbsp)[0]
    tau = kendalltau(logfc, logpbsp)[0]
    return {
        'name': name,
        'pwm': pwm_name,
        'batch': batch,
        'total': total,
        'positive': positive,
        'negative': negative,
        'auroc': auroc,
        'auprc': auprc,
        'r': r,
        'tau': tau,
        'metric': metric
    }


def parse_snpselex(name, batch):
    df = pd.DataFrame()
    try:
        df = pd.read_csv(BASE_PATH + f'/selex/batch{batch}/{name}.csv')
    except FileNotFoundError:
        return df
    df = df.apply(lambda row: parse_cols(row, batch), axis=1)
    return df


def parse_cols(row, batch):
    if batch == 1:
        chrom, pos = row['oligo'].split(':')
        row['#chr'] = chrom
        pos = tuple(map(int, pos.split('-')))
        assert pos[1] - pos[0] == 40
        row['pos'] = pos[0] + 20
        row['ID'] = row['rsid']
        return row[['#chr', 'pos', 'ID', 'ref', 'alt',
                    'pbs', 'pval', 'oligo_pval']]
    elif batch == 2:
        chrom, pos, ref, alt = row['snp'].split('_')
        row['#chr'] = chrom
        row['pos'] = int(pos)
        row['ref'] = ref
        row['alt'] = alt
        row['ID'] = row['snp']
        return row[['#chr', 'pos', 'ID', 'ref', 'alt',
                    'pbs', 'pval', 'oligo_pval']]


def add_ground_truth(row, batch):
    if row['pval'] < TP_THR and (row['oligo_pval'] < LH_THR or batch == 2):
        row['true'] = True
    elif row['pval'] > TN_THR and (row['oligo_pval'] < LH_THR or batch == 2):
        row['true'] = False
    else:
        row['true'] = None
    return row


def get_motif_length(pwm_path):
    with open(pwm_path, 'r') as f:
        return len(f.readlines()) - 1


def make_fasta_file(df, fasta, motif_length):
    if not os.path.isfile(fasta):
        if motif_length is None:
            print('No such PWM')
            sys.exit(0)
        with open(fasta, 'w') as outfile:
            df = df.apply(
                lambda row: extract_sarus(row, motif_length, outfile), axis=1)


def make_sarus_file(fasta, pvalues, pwm, thr):
    if not os.path.isfile(pvalues):
        output = get_pvalues(fasta, pwm, thr)
        with open(pvalues, 'wb') as out:
            out.write(output)


def extract_sarus(row, motif_length, outfile):
    chr = row['#chr']
    pos = row['pos']
    ref = row['ref']
    alt = row['alt']
    start = pos - motif_length
    end = pos + motif_length - 1
    id = f'{row["ID"]}_{row["pval"]:.5f}'
    coords = subprocess.Popen(['echo', f'{chr}\t{start}\t{end}'],
                              stdout=subprocess.PIPE)
    output = subprocess.check_output(['./bedtools.static', 'getfasta',
                                      '-fi', REF_PATH, '-bed', 'stdin'],
                                     stdin=coords.stdout)
    sequence = str(output, 'utf-8').split('\n')[1]
    if ref.upper() != sequence[motif_length - 1].upper():
        print(f'{id}: Failed')
    else:
        left = sequence[:motif_length-1].upper()
        right = sequence[motif_length:].upper()
        outfile.write(f'>{id}_ref\n{left}{ref}{right}\n')
        outfile.write(f'>{id}_alt\n{left}{alt}{right}\n')


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


def get_snp_dict(pvalues_path, motif_length):
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


def add_score_best_p(row, snp_dict, metric):
    ID = f'{row["ID"]}_{row["pval"]:.5f}'
    ref_p = snp_dict[ID]['ref']['p']
    alt_p = snp_dict[ID]['alt']['p']
    best_p = max(ref_p, alt_p)
    row['best_p'] = best_p
    row['logfc'] = alt_p - ref_p
    if row['pval'] != 0:
        row['logpbsp'] = -np.log10(row['pval']) * np.sign(-row['pbs'])
    else:
        row['logpbsp'] = None
    row['score'] = abs(row['logfc'])
    return row


if __name__ == '__main__':
    with open('selex_motifs.tsv', 'a') as results:
        tf_name = sys.argv[1].split('@')[0]
        pwm_name = '@'.join(sys.argv[1].split('@')[1:])
        batches = [1, 2]
        for batch in batches:
            result = main(tf_name, pwm_name, batch)
            row = [result['name'],
                   result['pwm'],
                   str(result['batch']),
                   str(result['total']),
                   str(result['positive']),
                   str(result['negative']),
                   str(result['auroc']),
                   str(result['auprc']),
                   str(result['r']),
                   str(result['tau'])]
            results.write('\t'.join(row) + '\n')
