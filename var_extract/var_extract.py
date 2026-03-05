#!/usr/bin/env python
# encoding: utf-8

# Reads BAM files and reports per-read variant calls

import pandas as pd
import pysam

from collections import defaultdict as ddict
import sys, os, glob, re


def print_log(msg, verbose=True):
    if verbose:
        print(msg)

def parse_vars(vars_fn, verbose=True):
    print_log(
            '>Reading variants input table', verbose)
    if not vars_fn.endswith('.csv'):
        sys.exit('Only csv format accepted as variants input table. Exiting')
    variants = pd.read_csv(vars_fn)
    if not {'pos', 'ref', 'alt'}.issubset(variants.columns):
        sys.exit((
            'At least one required column (pos, ref, alt) '
            'not present in variants input table. Exiting'))
        print_log(
                '>Variants input table read correctly',
                verbose)
    return variants

def filter_vars(variants, accept_indels=False, verbose=True):
    print_log(
            '>Filtering input variants', verbose)
    variants = variants.sort_values(
            by=['pos', 'ref', 'alt']).drop_duplicates()
    print_log(
            f'>{variants.shape[0]} distinct variants in input table',
            verbose)

    if not accept_indels:
        print_log(
                '>Ignoring INDELs', verbose)
        filtered_variants = variants[len(variants['ref']) == 1]
        filtered_variants = variants[len(variants['alt']) == 1]
    else:
        filtered_variants = variants

    if variants.shape[0] == filtered_variants.shape[0]:
        print_log('>All variants retained after filters', verbose)
    else:
        print_log(
            f'>{filtered_variants.shape[0]} retained after filters',
            verbose)
        print_log('>Filtered out variants:', verbose)
        print_log(f'{filtered_variants}', verbose)
    return filtered_variants

def read_alignment(als_fn, verbose=True):
    print_log(
        f'>Reading alignment file {als_fn}',
        verbose)
    als = pysam.AlignmentFile(als_fn, 'rb')
    print_log(
        '>Alignment file read successfully',
        verbose)
    return als

def extract_bases_pos(
        alignments, pos_ref,
        region = 'KY962518.1_looped_2120',
        loop_size = 2120,
        verbose = True):
    print_log(
        f'>Extracting sequence at variant position {pos_ref}',
        verbose)
    pstart = pos_ref + loop_size if pos_ref < 0 else pos_ref + loop_size - 1
    pend = pstart + 1
    bases = []
    for pcol in alignments.pileup(region, pstart, pend,
            ignore_orphans = False):
        if pcol.pos < pstart or pcol.pos >= pend:
            continue
        for pread in pcol.pileups:
            pname = pread.alignment.query_name
            indel = pread.indel
            pseq = pread.alignment.query_sequence
            if pread.is_del or pread.is_refskip:
                pbase = ''
                indel = 0
            else:
                ppos = pread.query_position
                if indel > 0:
                    ppos = pread.query_position
                    pbase = pseq[ppos:(ppos+indel+1)]
                else:
                    pbase = pseq[ppos]
            bases.append({
                'read_name' : pname,
                'pos_ref' : pos_ref,
                'obs_seq' : pbase,
                'indel_size' : indel
                })
    return pd.DataFrame(bases)

def extract_bases(alignments, variants,
        region = 'KY962518.1_looped_2120',
        loop_size = 2120,
        verbose = True):
    print_log(
        '>Extracting sequences at variant positions from alignment file',
        verbose)
    extracted_bases = []
    for pos_ref in set(variants['pos']):
        extracted_bases.append(
                extract_bases_pos(
                    alignments, pos_ref, region, loop_size, verbose))
    print_log(
        '>Sequences extracted', verbose)
    return pd.concat(extracted_bases)

def _parse_args():
    import argparse
    parser = argparse.ArgumentParser(
        description='Obtain per-read variant calls.')
    parser.add_argument('sample', metavar='S', type=str,
                               help='root name of sample to process')
    parser.add_argument('-r', '--region', type=str, default='KY962518.1_looped_2120',
                               help='contig name')
    parser.add_argument('-o', '--output', type=str,
                               help='path to output folder')
    parser.add_argument('-a', '--align', type=str,
                               help='path to alignments folder')
    parser.add_argument('-l', '--variants', type=str,
                               help='path to variants list in csv format')
    parser.add_argument('-v', '--verbose', action='store_true',
                               help='increase output verbosity')
    parser.add_argument('-i', '--indels', action='store_true',
                               help='include INDELs in the analysis')


    return parser.parse_args()

def main():

# 00.- Configuration

    args = _parse_args()
    root_fn = args.sample
    region = args.region
    if re.search('_looped_', region):
        loop_size = int(re.split('_', region)[-1].strip())
    else:
        loop_size = 0
    print(f"Starting execution for sample {root_fn}")

## Paths

    ofd = f"{args.output}/{root_fn}/"

    if not os.path.exists(ofd):
            os.makedirs(ofd)

    als_fd = f"{args.align}/{root_fn}/"
    als_fn = glob.glob(f"{als_fd}*.bam")[0]

    vars_fn = args.variants

## Options

    verbose = args.verbose
    indels = args.indels

# 01.- Read Variants

    variants = filter_vars(
            parse_vars(vars_fn, verbose = verbose),
            verbose = verbose, accept_indels = indels)

# 02.- Read alignment files

    als = read_alignment(als_fn, verbose = verbose)

# 03.- Extract bases from reads

    calls = pd.DataFrame(
            extract_bases(
                als,
                variants,
                region = region,
                loop_size = loop_size,
                verbose = verbose))

    calls = calls.sort_values(by=['read_name', 'pos_ref'])

    calls['read_num'] = calls.groupby('read_name').ngroup()

    print(f">Saving output for sample {root_fn}")
    reads = calls[['read_num', 'read_name']].drop_duplicates()
    reads.to_csv(f"{ofd}{root_fn}_reads.csv",
            index=False)

    calls = calls[['read_num', 'pos_ref', 'obs_seq', 'indel_size']]
    calls.to_csv(f"{ofd}{root_fn}_calls.csv",
            index=False)

if __name__ == "__main__":
    main()
