#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2022-07-24 17:05


import argparse
import collections
import os
import re
def parse_fasta(fasta_file):

    group_ctg_dict = collections.defaultdict(list)
    fa_dict = collections.defaultdict(list)    
    
    with open(fasta_file) as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith('>'):
                ctg = line.split()[0][1:]
                group = ctg.split('_')[0] +"_" + ctg.split('_')[1]
                assert ctg not in group_ctg_dict[group]
                group_ctg_dict[group].append(ctg)
            else:
                fa_dict[ctg].append(line.strip().upper())
    
    for ctg, seq_list in fa_dict.items():
        seq = ''.join(seq_list)
        # Using fixed 'GATC' here is ok. Becuase in the sorting step,
        # we don't acre about the restriction sites
        RE_sites = seq.count('GATC') #   AAGCTT
        fa_dict[ctg] = [seq, len(seq), RE_sites]

    return fa_dict, group_ctg_dict

def split_clm_file(clm_file, group_ctg_dict, ctg_group_dict, subdir):


    # make directory for clm splitting
    subdir = 'split_clms'
    os.mkdir(subdir)

    fp_dict = dict()

    for group in group_ctg_dict:
        fp_dict[group] = open('{}/{}.clm'.format(subdir, group), 'w')

    with open(clm_file) as f:
        for line in f:
            cols = line.split()
            ctg_1, ctg_2 = cols[0][:-1], cols[1][:-1]
            if ctg_1 in ctg_group_dict and ctg_2 in ctg_group_dict and ctg_group_dict[ctg_1] == ctg_group_dict[ctg_2]:
                fp_dict[ctg_group_dict[ctg_1]].write(line)

    for group, fp in fp_dict.items():
        fp.close()

def generate_group_files(fa_dict, group_ctg_dict):
    new_ctg_group_dict={}
    new_group_ctg_dict={}
    # sort by length
    for group, ctgs in group_ctg_dict.items():
        if group.lower().startswith("scaffold") or group.startswith("Contig") or group.startswith("tig") or re.match(r'^\d', group):
            continue
        with open('group_{}.txt'.format(group), 'w') as fout:
            fout.write('#Contig\tRECounts\tLength\n')
            sorted_ctgs = sorted(ctgs, key=lambda x: fa_dict[x][1], reverse=True)
            group_name='group_{}'.format(group)
            ctg_stat=[set(), 0]
            ctg_stat[0]= set(sorted_ctgs)
            for ctg in sorted_ctgs:
                new_ctg_group_dict[ctg] = group_name
                ctg_stat[1]+= fa_dict[ctg][1]
                fout.write('{}\t{}\t{}\n'.format(ctg, fa_dict[ctg][2], fa_dict[ctg][1]))
            new_group_ctg_dict[group_name] = ctg_stat
    return new_group_ctg_dict, new_ctg_group_dict

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('fasta', help='contig-level genome file in FASTA format')
    args = parser.parse_args()

    fa_dict, group_ctg_dict = parse_fasta(args.fasta)
    print(group_ctg_dict.keys())
    new_group_ctg_dict, new_ctg_group_dict=generate_group_files(fa_dict, group_ctg_dict)
    print(new_group_ctg_dict.keys())
    print(new_ctg_group_dict.keys())
    split_clm_file("/home/chenzh/HapHiC/data/XinJiangDaYe/with_parameter_in_article/01.cluster/paired_links.clm", new_group_ctg_dict, new_ctg_group_dict, 'hc_groups')

if __name__ == '__main__':
    main()

