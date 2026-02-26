import sys
import os
import pickle
import argparse
import logging
import time
import pysam
import gzip

from collections import defaultdict
from sklearn.cluster import AgglomerativeClustering

import warnings
warnings.filterwarnings('ignore', category=FutureWarning)

from ConHiC_cluster import (
    parse_fasta, parse_gfa, stat_fragments,
    remove_allelic_HiC_links, dict_to_matrix
)

from _version import __version__, __update_time__

logging.basicConfig(
    format='%(asctime)s <%(filename)s> [%(funcName)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

pysam.set_verbosity(0)


# ─────────────────────────────────────────────────────────────
# 输入文件解析
# ─────────────────────────────────────────────────────────────

def parse_pickle(fa_dict, pickle_file):
    """
    从 pickle 文件中读取 contig 间完整链接字典，并生成辅助数据结构。

    参数
    ----
    fa_dict : dict
        parse_fasta 返回的 contig 信息字典
    pickle_file : str
        full_links.pkl 文件路径

    返回
    ----
    tuple
        (full_link_dict, sorted_ctg_list, RE_site_dict)
        其中 sorted_ctg_list 按 contig 长度降序排列
    """
    logger.info('读取 pickle 文件…')

    with open(pickle_file, 'rb') as fh:
        full_link_dict = pickle.load(fh)

    sorted_ctg_list = sorted(
        [(ctg, fa_dict[ctg][1]) for ctg in fa_dict],
        key=lambda x: x[1], reverse=True
    )
    RE_site_dict = {ctg: info[2] for ctg, info in fa_dict.items()}

    return full_link_dict, sorted_ctg_list, RE_site_dict

def parse_clusters(clusters_file, RE_site_dict, fa_dict, min_group_len):
    """
    解析 MCL 聚类步骤输出的 .clusters.txt 文件，建立 contig→组 的映射关系。

    长度低于 min_group_len 的组将被忽略，其中的 contig 随后会被重新分配。

    参数
    ----
    clusters_file : str
        聚类文件路径
    RE_site_dict : dict
        {contig: RE 位点数}
    fa_dict : dict
        contig 信息字典
    min_group_len : float
        组的最小长度（Mbp），低于此值的组将被丢弃

    返回
    ----
    tuple
        (ctg_group_dict, group_RE_dict)
    """
    logger.info('解析 .clusters.txt 文件…')

    ctg_group_dict = {}
    group_RE_dict = {}

    with open(clusters_file) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            cols = line.split()
            group = cols[0]
            if min_group_len and sum(fa_dict[c][1] for c in cols[2:]) / 1_000_000 < min_group_len:
                continue
            group_RE_dict[group] = 1
            for ctg in cols[2:]:
                ctg_group_dict[ctg] = group
                group_RE_dict[group] += RE_site_dict[ctg] - 1

    return ctg_group_dict, group_RE_dict


def parse_assembly(assembly_file, RE_site_dict, fa_dict, min_group_len):
    """
    解析 .assembly 格式文件，建立 contig→组 的映射关系。

    参数
    ----
    assembly_file : str
        .assembly 文件路径
    RE_site_dict : dict
        {contig: RE 位点数}
    fa_dict : dict
        contig 信息字典
    min_group_len : float
        组的最小长度（Mbp）

    返回
    ----
    tuple
        (ctg_group_dict, group_RE_dict)
    """
    logger.info('解析 .assembly 文件…')

    ctg_dict = {}
    ctg_group_dict = {}
    group_RE_dict = {}
    n = 0

    with open(assembly_file) as fh:
        for line in fh:
            if not line.strip():
                continue
            cols = line.split()
            if line.startswith('>'):
                ctg_dict[cols[1]] = cols[0][1:]
            else:
                n += 1
                group = 'group{}'.format(n)
                if min_group_len and sum(fa_dict[ctg_dict[num.strip('-')]][1] for num in cols) / 1_000_000 < min_group_len:
                    continue
                group_RE_dict[group] = 1
                for num in cols:
                    ctg = ctg_dict[num.strip('-')]
                    ctg_group_dict[ctg] = group
                    group_RE_dict[group] += RE_site_dict[ctg] - 1

    return ctg_group_dict, group_RE_dict


# ─────────────────────────────────────────────────────────────
# 链接解析与统计
# ─────────────────────────────────────────────────────────────

def add_ungrouped_ctgs(fa_dict, ctg_group_dict):
    """
    将 fa_dict 中所有尚未分组的 contig 标记为 'ungrouped'，
    并返回已分组 contig 的集合（用于后续层次聚类步骤）。

    参数
    ----
    fa_dict : dict
        contig 信息字典
    ctg_group_dict : dict
        {contig: 组名}（将被原地更新）

    返回
    ----
    set
        已分组 contig 的集合
    """
    grouped = set()
    for ctg in fa_dict:
        if ctg not in ctg_group_dict:
            ctg_group_dict[ctg] = 'ungrouped'
        else:
            grouped.add(ctg)
    return grouped


def parse_link_dict(link_dict, ctg_group_dict, normalize_by_nlinks=False):
    """
    统计每条 contig 与各组之间的 Hi-C 链接总数，并可选地按总链接数进行归一化。

    参数
    ----
    link_dict : dict
        {(contig_i, contig_j): 链接数}
    ctg_group_dict : dict
        {contig: 组名}
    normalize_by_nlinks : bool
        是否按各 contig 总链接数的几何平均值进行归一化

    返回
    ----
    tuple
        (ctg_group_link_dict, linked_ctg_dict)
        ctg_group_link_dict : {contig: {组名: 链接数}}
        linked_ctg_dict : {contig: 与其有链接的 contig 集合}
    """
    ctg_group_link_dict = defaultdict(dict)
    linked_ctg_dict = defaultdict(set)
    ctg_total_links = defaultdict(int)

    def register(ctg, group, links):
        if group != 'ungrouped':
            ctg_group_link_dict[ctg][group] = ctg_group_link_dict[ctg].get(group, 0) + links

    if normalize_by_nlinks:
        total_raw = 0
        for (ci, cj), lk in link_dict.items():
            ctg_total_links[ci] += lk
            ctg_total_links[cj] += lk
            total_raw += lk

    total_norm = 0
    for (ci, cj), lk in link_dict.items():
        if normalize_by_nlinks:
            lk /= (ctg_total_links[ci] * ctg_total_links[cj]) ** 0.5
            link_dict[(ci, cj)] = lk
            total_norm += lk
        else:
            gi, gj = ctg_group_dict[ci], ctg_group_dict[cj]
            register(ci, gj, lk)
            register(cj, gi, lk)
            linked_ctg_dict[ci].add(cj)
            linked_ctg_dict[cj].add(ci)

    if normalize_by_nlinks:
        scale = total_raw / total_norm
        for (ci, cj), lk in link_dict.items():
            lk *= scale
            link_dict[(ci, cj)] = lk
            gi, gj = ctg_group_dict[ci], ctg_group_dict[cj]
            register(ci, gj, lk)
            register(cj, gi, lk)
            linked_ctg_dict[ci].add(cj)
            linked_ctg_dict[cj].add(ci)

    return ctg_group_link_dict, linked_ctg_dict


# ─────────────────────────────────────────────────────────────
# 重分配核心逻辑
# ─────────────────────────────────────────────────────────────

def run_reassignment(
        sorted_ctg_list, ctg_group_link_dict, ctg_group_dict,
        full_link_dict, linked_ctg_dict, fa_dict, RE_site_dict,
        group_RE_dict, max_ctg_len, min_RE_sites, min_links,
        min_link_density, min_density_ratio, ambiguous_cutoff,
        min_group_len, whitelist, nround):
    """
    对所有 contig 依长度从大到小遍历，依据 Hi-C 链接密度将其分配或保留至最优组。

    本函数支持两种模式：
    - nround > 0：正式重分配轮次，同时允许将已分组的 contig 迁移至更优的组；
    - nround == 0：额外救援轮次，仅将未分组 contig 纳入最优组，不进行组间迁移。

    过滤条件（任意一条成立则跳过当前 contig）：
    （1）RE 位点数低于阈值；
    （2）与任意组的链接数均低于阈值；
    （3）正式轮次中，次优组链接比例超过模糊截断值；
    （4）链接密度低于阈值。

    参数
    ----
    sorted_ctg_list : list
        按长度降序排列的 (contig, 长度) 列表
    ctg_group_link_dict : dict
        {contig: {组名: 链接数}}，将被原地更新
    ctg_group_dict : dict
        {contig: 当前所属组名}，将被原地更新
    full_link_dict : dict
        完整的 contig 间链接字典
    linked_ctg_dict : dict
        {contig: 与其有链接的 contig 集合}
    fa_dict : dict
        contig 信息字典
    RE_site_dict : dict
        {contig: RE 位点数}
    group_RE_dict : dict
        {组名: 组内 RE 位点总数}，将被原地更新
    max_ctg_len : float
        允许重分配的最大 contig 长度（kbp）
    min_RE_sites : int
        RE 位点数下限
    min_links : int
        最小链接数
    min_link_density : float
        最小链接密度
    min_density_ratio : float
        最小链接密度比值（最优组密度 / 其他组平均密度）
    ambiguous_cutoff : float
        模糊 contig 判定阈值（次优 / 最优链接比）
    min_group_len : float
        组的最小长度（Mbp）；小组将被解散
    whitelist : set
        白名单 contig 集合，豁免过滤
    nround : int
        当前轮次编号（0 表示额外救援轮次）
    """

    def transfer(ctg, target_group):
        """将 ctg 从当前组迁移至 target_group，并同步更新链接统计。"""
        src_group = ctg_group_dict[ctg]
        ctg_group_dict[ctg] = target_group
        for nbr in linked_ctg_dict[ctg]:
            gl = ctg_group_link_dict[nbr]
            lk = full_link_dict[tuple(sorted([ctg, nbr]))]
            if src_group != 'ungrouped':
                gl[src_group] = gl.get(src_group, 0) - lk
            if target_group != 'ungrouped':
                gl[target_group] = gl.get(target_group, 0) + lk

    def link_density(ctg, group, current_group, links):
        """计算 contig 与目标组之间的链接密度。"""
        re_g = group_RE_dict[group]
        if group == current_group:
            return links / re_g
        return links / (re_g + RE_site_dict[ctg] - 1)

    round_name = 'round{}'.format(nround) if nround else 'additional_rescue'
    if nround:
        logger.info('执行重分配（{}）…'.format(round_name))
    else:
        logger.info('执行额外救援轮次…')

    # 解散过小的组（仅从第 2 轮起生效）
    if min_group_len and nround > 1:
        grp_len = defaultdict(int)
        for ctg, grp in ctg_group_dict.items():
            if grp != 'ungrouped':
                grp_len[grp] += fa_dict[ctg][1]

        dismissed = {g for g, l in grp_len.items() if l / 1_000_000 < min_group_len}
        logger.info('[reassignment::{}] 已解散的组：{}'.format(round_name, dismissed))
        for ctg in ctg_group_dict:
            if ctg_group_dict[ctg] in dismissed:
                ctg_group_dict[ctg] = 'ungrouped'
            for dg in dismissed:
                ctg_group_link_dict[ctg][dg] = 0

    result = defaultdict(int)
    skip_count = 0
    filter_counts = [0, 0, 0, 0, 0]  # [RE, links, ambiguous, density, ratio]

    with open('output.txt', 'w') as fout:
        for ctg, ctg_len in sorted_ctg_list:
            src_group = ctg_group_dict[ctg]
            grp_links = ctg_group_link_dict[ctg]
            filtered = False

            # 条件 1：RE 位点过滤
            if (RE_site_dict[ctg] - 1 < min_RE_sites and ctg not in whitelist) or not grp_links:
                filter_counts[0] += 1
                max_group, max_links = 'ungrouped', 0
                logger.debug('[reassignment::{}] {} 未被救援（RE 位点不足，当前组：{}）'.format(
                    round_name, ctg, src_group))
                result['not_rescued'] += 1
                filtered = True
            else:
                sorted_gl = sorted(grp_links.items(), key=lambda x: x[1], reverse=True)
                max_group, max_links = sorted_gl[0]
                second_links = sorted_gl[1][1] if len(sorted_gl) > 1 else 0

                # 条件 2：最优链接数过滤
                if max_links < min_links and ctg not in whitelist:
                    filter_counts[1] += 1
                    max_group, max_links = 'ungrouped', 0
                    logger.debug('[reassignment::{}] {} 未被救援（链接数不足，当前组：{}）'.format(
                        round_name, ctg, src_group))
                    result['not_rescued'] += 1
                    filtered = True

                # 条件 3：模糊 contig 过滤（仅正式轮次）
                elif nround and second_links / max_links >= ambiguous_cutoff and ctg not in whitelist:
                    filter_counts[2] += 1
                    max_group, max_links = 'ungrouped', 0
                    logger.debug('[reassignment::{}] {} 未被救援（链接模糊，当前组：{}）'.format(
                        round_name, ctg, src_group))
                    result['not_rescued'] += 1
                    filtered = True
                else:
                    max_ld = link_density(ctg, max_group, src_group, max_links)

                    # 条件 4：最大链接密度过滤
                    if max_ld < min_link_density and ctg not in whitelist:
                        filter_counts[3] += 1
                        max_group, max_links = 'ungrouped', 0
                        logger.debug('[reassignment::{}] {} 未被救援（链接密度过低，当前组：{}）'.format(
                            round_name, ctg, src_group))
                        result['not_rescued'] += 1
                        filtered = True
                    else:
                        # 计算其他所有组的平均链接密度
                        other_ld_sum = sum(
                            link_density(ctg, g, src_group, lk)
                            for g, lk in sorted_gl[1:]
                        )
                        n_others = len(group_RE_dict) - 1
                        avg_other_ld = other_ld_sum / n_others if other_ld_sum else 1_000_000_000

            if filtered:
                skip_count += 1
                continue

            if max_group == 'ungrouped':
                assert False
                logger.debug('[reassignment::{}] {} 保持不变（{}）'.format(round_name, ctg, src_group))
                result['consistent'] += 1
            else:
                ld_ratio = max_ld / avg_other_ld

                # 未分组 contig 的救援
                if src_group == 'ungrouped':
                    if ld_ratio >= min_density_ratio:
                        transfer(ctg, max_group)
                        group_RE_dict[max_group] += RE_site_dict[ctg] - 1
                        logger.debug('[reassignment::{}] {} 已救援（ungrouped → {}）'.format(
                            round_name, ctg, max_group))
                        result['rescued'] += 1
                    else:
                        filter_counts[4] += 1
                        logger.debug('[reassignment::{}] {} 未被救援（密度比值过低）'.format(round_name, ctg))
                        result['not_rescued'] += 1

                # 当前组即为最优组，保持不变
                elif src_group in grp_links and grp_links[src_group] == max_links:
                    logger.debug('[reassignment::{}] {} 保持不变（{}）'.format(round_name, ctg, src_group))
                    result['consistent'] += 1

                # 正式轮次中的组间重分配
                elif nround and ctg_len <= max_ctg_len * 1000 and ld_ratio >= min_density_ratio:
                    transfer(ctg, max_group)
                    group_RE_dict[src_group] -= RE_site_dict[ctg] - 1
                    group_RE_dict[max_group] += RE_site_dict[ctg] - 1
                    logger.debug('[reassignment::{}] {} 已重分配（{} → {}）'.format(
                        round_name, ctg, src_group, max_group))
                    result['reassigned'] += 1
                else:
                    logger.debug('[reassignment::{}] {} 保持不变（{}）'.format(round_name, ctg, src_group))
                    result['consistent'] += 1

        print('跳过总数：', skip_count)
        print('各过滤条件计数：', filter_counts)
        logger.info(
            '[result::{}] 总计：{}，保持不变：{}，救援：{}，重分配：{}，未救援：{}'.format(
                round_name, len(sorted_ctg_list),
                result['consistent'], result['rescued'],
                result['reassigned'], result['not_rescued']
            )
        )


# ─────────────────────────────────────────────────────────────
# 聚类结果统计与输出
# ─────────────────────────────────────────────────────────────

def stat_clusters(ctg_group_dict, fa_dict, grouped_ctgs):
    """
    统计重分配后每个组内的 contig 集合、总长度及高置信 RE 位点数。

    高置信 contig 指在 MCL 聚类步骤中已被分组的 contig，
    其 RE 位点数将用于后续层次聚类步骤的密度计算。

    参数
    ----
    ctg_group_dict : dict
        {contig: 组名}
    fa_dict : dict
        contig 信息字典
    grouped_ctgs : set
        MCL 聚类阶段已分组的高置信 contig 集合

    返回
    ----
    tuple
        (group_ctg_dict, group_hiconf_RE_dict)
        group_ctg_dict : {组名: [{(contig, 长度)}, 总长度]}
        group_hiconf_RE_dict : {组名: 高置信 RE 位点总数}
    """
    group_ctg_dict = {}
    group_hiconf_RE_dict = {}

    for ctg, group in ctg_group_dict.items():
        if group == 'ungrouped':
            continue
        clen = fa_dict[ctg][1]
        re_contrib = fa_dict[ctg][2] - 1 if ctg in grouped_ctgs else 0

        if group in group_ctg_dict:
            group_ctg_dict[group][0].add((ctg, clen))
            group_ctg_dict[group][1] += clen
            group_hiconf_RE_dict[group] += re_contrib
        else:
            group_ctg_dict[group] = [{(ctg, clen)}, clen]
            group_hiconf_RE_dict[group] = re_contrib

    return group_ctg_dict, group_hiconf_RE_dict


def clusters_output(group_ctg_dict, fa_dict, out_prefix):
    """
    将重分配或层次聚类结果写入分组文件，并返回更新后的映射字典。

    每个组生成一个独立的 TXT 文件，同时生成汇总的 _clusters.txt 文件。

    参数
    ----
    group_ctg_dict : dict
        {组名: [{(contig, 长度)}, 总长度]}
    fa_dict : dict
        contig 信息字典
    out_prefix : str
        输出前缀，'reassigned' 或 'hc'

    返回
    ----
    tuple
        (new_ctg_group_dict, new_group_ctg_dict, new_old_group_dict)
    """
    new_ctg_group_dict = {}
    new_old_group_dict = {}
    new_group_ctg_dict = {}

    sorted_clusters = sorted(group_ctg_dict.items(), key=lambda x: x[1][1], reverse=True)
    subdir = 'reassigned_groups' if out_prefix == 'reassigned' else 'hc_groups'

    with open('{}/{}_clusters.txt'.format(subdir, out_prefix), 'w') as fclusters:
        fclusters.write('#Group\tnContigs\tContigs\n')
        for n, (group, stat) in enumerate(sorted_clusters, 1):
            gname = 'group{}_{}bp'.format(n, stat[1])
            new_old_group_dict[gname] = group
            new_group_ctg_dict[gname] = stat

            sorted_ctgs = [c for c, _ in sorted(stat[0], key=lambda x: x[1], reverse=True)]
            fclusters.write('{}\t{}\t{}\n'.format(gname, len(stat[0]), ' '.join(sorted_ctgs)))

            with open('{}/{}_{}.txt'.format(subdir, out_prefix, gname), 'w') as fgrp:
                fgrp.write('#Contig\tRECounts\tLength\n')
                for ctg in sorted_ctgs:
                    fgrp.write('{}\t{}\t{}\n'.format(ctg, fa_dict[ctg][2], fa_dict[ctg][1]))
                    new_ctg_group_dict[ctg] = gname

    return new_ctg_group_dict, new_group_ctg_dict, new_old_group_dict


# ─────────────────────────────────────────────────────────────
# 层次聚类（合并小组）
# ─────────────────────────────────────────────────────────────

def agglomerative_hierarchical_clustering(
        full_link_dict, grouped_ctgs, new_ctg_group_dict,
        group_hiconf_RE_dict, new_old_group_dict, nclusters,
        normalize_by_nlinks=False):
    """
    对重分配后的各组执行凝聚层次聚类，将组数压缩至预期染色体数。

    以组间链接密度的倒置作为距离度量，使用平均链接策略合并最相近的组。

    参数
    ----
    full_link_dict : dict
        完整的 contig 间链接字典
    grouped_ctgs : set
        高置信 contig 集合
    new_ctg_group_dict : dict
        {contig: 重分配后的组名}
    group_hiconf_RE_dict : dict
        {组名: 高置信 RE 位点总数}
    new_old_group_dict : dict
        {新组名: 原组名}
    nclusters : int
        目标聚类数（通常等于预期染色体数）
    normalize_by_nlinks : bool
        是否按总链接数归一化

    返回
    ----
    dict
        {新聚类编号: [旧组名列表]}
    """
    logger.info('执行凝聚层次聚类（合并小组）…')

    group_link_dict = defaultdict(int)

    for (ci, cj), lk in full_link_dict.items():
        if ci not in grouped_ctgs or cj not in grouped_ctgs:
            continue
        if ci not in new_ctg_group_dict or cj not in new_ctg_group_dict:
            continue
        gi, gj = new_ctg_group_dict[ci], new_ctg_group_dict[cj]
        if gi == gj:
            continue
        group_link_dict[tuple(sorted([gi, gj]))] += lk

    # 计算组间链接密度
    with open('hc_groups/group_group_links.txt', 'w') as fout:
        fout.write('group1\tgroup2\tlinks\tlink_density\n')
        max_ld = 0

        if normalize_by_nlinks:
            grp_total = defaultdict(int)
            for (gi, gj), lk in group_link_dict.items():
                grp_total[gi] += lk
                grp_total[gj] += lk

        for (gi, gj), lk in group_link_dict.items():
            re_i = group_hiconf_RE_dict[new_old_group_dict[gi]]
            re_j = group_hiconf_RE_dict[new_old_group_dict[gj]]
            if normalize_by_nlinks:
                ld = lk / (grp_total[gi] * grp_total[gj])
            else:
                ld = lk / (re_i * re_j)
            max_ld = max(max_ld, ld)
            group_link_dict[(gi, gj)] = ld
            fout.write('{}\t{}\t{}\t{}\n'.format(gi, gj, lk, ld))

    # 链接密度 → 距离矩阵
    mat, grp_idx = dict_to_matrix(group_link_dict, set(new_old_group_dict.keys()))
    idx_grp = {i: g for g, i in grp_idx.items()}
    dist_mat = max_ld - mat

    # 凝聚层次聚类
    if 'affinity' in AgglomerativeClustering._get_param_names():
        clust = AgglomerativeClustering(
            n_clusters=nclusters, affinity='precomputed',
            linkage='average', distance_threshold=None
        )
    else:
        logger.info('scikit-learn 版本较新，使用 metric 参数代替 affinity')
        assert 'metric' in AgglomerativeClustering._get_param_names()
        clust = AgglomerativeClustering(
            n_clusters=nclusters, metric='precomputed',
            linkage='average', distance_threshold=None
        )

    labels = clust.fit_predict(dist_mat)

    hc_dict = defaultdict(list)
    for idx, label in enumerate(labels):
        hc_dict[label].append(idx_grp[idx])

    return hc_dict


def stat_hc_clusters(group_ctg_dict, hc_cluster_dict):
    """
    合并各组的 contig 信息，生成层次聚类后的组统计字典，并输出结果文件。

    参数
    ----
    group_ctg_dict : dict
        重分配阶段的 {组名: [{(contig, 长度)}, 总长度]}
    hc_cluster_dict : dict
        {层次聚类编号: [旧组名列表]}

    返回
    ----
    dict
        {层次聚类编号: [{(contig, 长度)}, 总长度]}
    """
    hc_group_dict = {}

    with open('hc_groups/hc_result.txt', 'w') as fout:
        fout.write('hc_id\treassigned_groups\n')
        for new_id, old_groups in hc_cluster_dict.items():
            fout.write('{}\t{}\n'.format(new_id, ' '.join(old_groups)))
            hc_group_dict[new_id] = [set(), 0]
            for og in old_groups:
                hc_group_dict[new_id][0] |= group_ctg_dict[og][0]
                hc_group_dict[new_id][1] += group_ctg_dict[og][1]

    return hc_group_dict


# ─────────────────────────────────────────────────────────────
# CLM 文件分割与快速视图
# ─────────────────────────────────────────────────────────────

def split_clm_file(clm_file, group_ctg_dict, ctg_group_dict, subdir):
    """
    将全局 CLM 文件按分组拆分为各组独立的子文件，供后续排序步骤使用。

    同时生成最终分组目录（final_groups），并以符号链接指向对应的分组文件。

    参数
    ----
    clm_file : str
        paired_links.clm 文件路径
    group_ctg_dict : dict
        {组名: contig 信息}
    ctg_group_dict : dict
        {contig: 组名}
    subdir : str
        来源子目录，'reassigned_groups' 或 'hc_groups'
    """
    logger.info('按组拆分 CLM 文件…')

    final_dir = 'final_groups'
    os.mkdir(final_dir)

    prefix = 'reassigned' if subdir == 'reassigned_groups' else 'hc'

    for group in group_ctg_dict:
        os.symlink(
            '../{0}/{1}_{2}.txt'.format(subdir, prefix, group),
            '{0}/{1}.txt'.format(final_dir, group)
        )
    os.symlink(
        '../{0}/{1}_clusters.txt'.format(subdir, prefix),
        '{0}/final_clusters.txt'.format(final_dir)
    )

    split_dir = 'split_clms'
    os.mkdir(split_dir)

    fps = {grp: open('{}/{}.clm'.format(split_dir, grp), 'w') for grp in group_ctg_dict}

    with open(clm_file) as fh:
        for line in fh:
            cols = line.split()
            c1, c2 = cols[0][:-1], cols[1][:-1]
            if (c1 in ctg_group_dict and c2 in ctg_group_dict
                    and ctg_group_dict[c1] == ctg_group_dict[c2]):
                fps[ctg_group_dict[c1]].write(line)

    for fp in fps.values():
        fp.close()


def mock_clusters_file(fa_dicts, total_lens, final_dir):
    """为快速视图模式生成占位 clusters 文件。"""
    with open('{}/final_clusters.txt'.format(final_dir), 'w') as fout:
        fout.write('#Group\tnContigs\tContigs\n')
        for n, (fa_dict, tlen) in enumerate(zip(fa_dicts, total_lens), 1):
            fout.write('group{}_{}bp\t{}\t{}\n'.format(
                n, tlen, len(fa_dict), ' '.join(fa_dict.keys())
            ))


def mock_group_file(fa_dicts, total_lens, final_dir):
    """为快速视图模式生成各组占位 TXT 文件。"""
    for n, (fa_dict, tlen) in enumerate(zip(fa_dicts, total_lens), 1):
        with open('{}/group{}_{}bp.txt'.format(final_dir, n, tlen), 'w') as fout:
            fout.write('#Contig\tRECounts\tLength\n')
            for ctg, info in fa_dict.items():
                fout.write('{}\t{}\t{}\n'.format(ctg, info[2], info[1]))


# ─────────────────────────────────────────────────────────────
# 参数解析
# ─────────────────────────────────────────────────────────────

def parse_arguments():
    """解析命令行参数并返回 Namespace 对象。"""

    parser = argparse.ArgumentParser(prog='conhic reassign')

    # 输入文件与流程控制
    inp = parser.add_argument_group('>>> 输入文件与流程控制参数')
    inp.add_argument('fasta',    help='草图基因组 FASTA 文件')
    inp.add_argument('links',    help='聚类步骤生成的 full_links.pkl，或 BAM/pairs 格式的 Hi-C 比对文件（请勿按坐标排序）')
    inp.add_argument('clusters', help='*.clusters.txt 或 *.assembly 聚类结果文件')
    inp.add_argument('clm',      help='聚类步骤生成的 paired_links.clm；重分配后将按组拆分以供排序步骤使用')
    inp.add_argument(
        '--RE', default='GATC',
        help='限制酶识别位点，多位点用逗号分隔，默认：%(default)s')
    inp.add_argument(
        '--quick_view', default=False, action='store_true',
        help='快速视图模式：跳过聚类与重分配，直接对全部 contig 进行快速排序，默认：%(default)s')
    inp.add_argument(
        '--gfa', default=None,
        help='（实验性）phased hifiasm 的 GFA 文件，多文件以逗号分隔，默认：%(default)s')

    # 重分配与救援参数
    rea = parser.add_argument_group('>>> 重分配与救援参数')
    rea.add_argument(
        '--min_group_len', type=float, default=5,
        help='组的最小长度（Mbp），低于此值的组将被解散，其中 contig 将被重新分配，默认：%(default)s')
    rea.add_argument(
        '--max_ctg_len', type=float, default=10000,
        help='允许重分配的最大 contig 长度（kbp），超过此长度的 contig 仅救援不重分配，默认：%(default)s')
    rea.add_argument(
        '--min_RE_sites', type=int, default=25,
        help='最小 RE 位点数，低于此值的 contig 将被移至 ungrouped，默认：%(default)s')
    rea.add_argument(
        '--min_links', type=int, default=25,
        help='最小 Hi-C 链接数，默认：%(default)s')
    rea.add_argument(
        '--min_link_density', type=float, default=0.0001,
        help='最小 Hi-C 链接密度（链接数 / RE 位点数），默认：%(default)s')
    rea.add_argument(
        '--min_density_ratio', type=float, default=4,
        help='最小链接密度比值（最优组密度 / 其他组平均密度），默认：%(default)s')
    rea.add_argument(
        '--ambiguous_cutoff', type=float, default=0.6,
        help='模糊 contig 判定阈值（次优 / 最优链接比），达到此比值的 contig 不参与重分配，默认：%(default)s')
    rea.add_argument(
        '--reassign_nrounds', type=int, default=0,
        help='最大重分配轮数，结果收敛时提前终止，默认：%(default)s')
    rea.add_argument(
        '--normalize_by_nlinks', default=False, action='store_true',
        help='按总链接数对 contig 间和组间链接进行归一化，默认：%(default)s')
    rea.add_argument(
        '--nclusters', type=int, default=0,
        help='重分配后执行凝聚层次聚类以合并小组，值通常设为预期染色体数；0 表示禁用，默认：%(default)s')
    rea.add_argument(
        '--no_additional_rescue', default=True, action='store_true',
        help='跳过额外救援轮次，默认：%(default)s')

    # 链接过滤参数（仅适用于 BAM/pairs 输入）
    flt = parser.add_argument_group('>>> Hi-C 链接过滤参数（仅对 BAM/pairs 格式输入生效）')
    flt.add_argument(
        '--remove_allelic_links', type=int, default=0,
        help='识别并移除等位 contig 之间的 Hi-C 链接，值为倍性数（≥2），默认禁用')
    flt.add_argument(
        '--concordance_ratio_cutoff', type=float, default=0.2,
        help='等位 contig 识别的一致性比率阈值，默认：%(default)s')
    flt.add_argument(
        '--nwindows', type=int, default=50,
        help='一致性比率计算中的分窗数目，默认：%(default)s')
    flt.add_argument(
        '--max_read_pairs', type=int, default=200,
        help='等位 contig 识别所用的最大 read 对数，默认：%(default)s')
    flt.add_argument(
        '--min_read_pairs', type=int, default=20,
        help='等位 contig 识别所需的最小 read 对数，默认：%(default)s')

    # 性能参数
    perf = parser.add_argument_group('>>> 性能参数')
    perf.add_argument(
        '--threads', type=int, default=8,
        help='读取 BAM 文件的线程数，默认：%(default)s')

    # 日志参数
    log = parser.add_argument_group('>>> 日志参数')
    log.add_argument(
        '--verbose', default=False, action='store_true',
        help='输出详细日志，默认：%(default)s')

    args = parser.parse_args()

    if not (args.links.endswith('.bam') or args.links.endswith('.pkl')
            or args.links.endswith('.pairs') or args.links.endswith('.pairs.gz')):
        logger.error('links 参数应以 .bam、.pkl、.pairs 或 .pairs.gz 结尾')
        raise RuntimeError('参数校验失败')

    if not args.clusters.endswith('.clusters.txt') and not args.clusters.endswith('.assembly'):
        logger.error('clusters 参数应以 .clusters.txt 或 .assembly 结尾')
        raise RuntimeError('参数校验失败')

    return args


# ─────────────────────────────────────────────────────────────
# 主流程
# ─────────────────────────────────────────────────────────────

def run(args, log_file=None, random_times=10):
    """
    ConHiC 重分配主流程。

    依次执行：读取 FASTA → 解析链接 → 解析聚类结果 → 多轮重分配与救援
    → （可选）凝聚层次聚类 → CLM 文件拆分。

    在 random_times > 1 时，将对多组聚类结果分别执行重分配，
    每组结果保存于以编号命名的子目录中。

    参数
    ----
    args : argparse.Namespace
        命令行参数
    log_file : str, optional
        附加日志文件路径
    random_times : int
        并行处理的聚类结果组数
    """
    if log_file:
        fh = logging.FileHandler(log_file, 'w')
        fh.setFormatter(logging.Formatter(
            fmt='%(asctime)s <%(filename)s> [%(funcName)s] %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        ))
        logger.addHandler(fh)

    start_time = time.time()
    logger.info('ConHiC 启动，版本：{} （更新日期：{}）'.format(__version__, __update_time__))
    logger.info('Python 版本：{}'.format(sys.version.replace('\n', '')))
    logger.info('命令：{}'.format(' '.join(sys.argv)))

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    # 读取基因组 FASTA
    fa_dict = parse_fasta(args.fasta, RE=args.RE, logger=logger)

    # 快速视图模式
    if args.quick_view:
        if args.gfa:
            gfa_list = args.gfa.split(',')
            read_depth_dict = parse_gfa(gfa_list, fa_dict)
        else:
            gfa_list = []
            read_depth_dict = {}

        final_dir = 'final_groups'
        os.mkdir(final_dir)

        if len(gfa_list) <= 1:
            total_len = sum(info[1] for info in fa_dict.values())
            mock_clusters_file((fa_dict,), (total_len,), final_dir)
            mock_group_file((fa_dict,), (total_len,), final_dir)
        else:
            hap_ctg = defaultdict(set)
            for ctg, (hap, _) in read_depth_dict.items():
                if ctg in fa_dict:
                    hap_ctg[hap].add(ctg)
            fa_dicts, total_lens = [], []
            for hap, ctgs in hap_ctg.items():
                hd = {c: fa_dict[c] for c in ctgs}
                fa_dicts.append(hd)
                total_lens.append(sum(info[1] for info in hd.values()))
            mock_clusters_file(fa_dicts, total_lens, final_dir)
            mock_group_file(fa_dicts, total_lens, final_dir)

        logger.info('程序完成，耗时 {}s'.format(time.time() - start_time))
        return None

    # 解析链接信息
    if args.links.endswith('.pkl'):
        if args.remove_allelic_links:
            logger.info('输入为 pickle 文件，--remove_allelic_links 在重分配步骤中不生效')
        full_link_dict, sorted_ctg_list, RE_site_dict = parse_pickle(fa_dict, args.links)

    # 对每组聚类结果分别执行重分配
    for i in range(random_times):
        os.makedirs(str(i), exist_ok=True)
        os.chdir(str(i))

        # 解析聚类文件
        cluster_path = args.clusters[i] if isinstance(args.clusters, list) else args.clusters
        if cluster_path.endswith('.clusters.txt'):
            ctg_group_dict, group_RE_dict = parse_clusters(
                cluster_path, RE_site_dict, fa_dict, args.min_group_len)
        else:
            ctg_group_dict, group_RE_dict = parse_assembly(
                cluster_path, RE_site_dict, fa_dict, args.min_group_len)

        # 标记未分组 contig
        grouped_ctgs = add_ungrouped_ctgs(fa_dict, ctg_group_dict)

        # 构建 contig → 组 的链接统计
        ctg_group_link_dict, linked_ctg_dict = parse_link_dict(
            full_link_dict, ctg_group_dict,
            normalize_by_nlinks=args.normalize_by_nlinks
        )

        logger.info('文件解析与数据准备完成，耗时 {}s'.format(time.time() - start_time))

        whitelist = getattr(args, 'whitelist', set())

        # 多轮重分配迭代
        for n in range(args.reassign_nrounds):
            run_reassignment(
                sorted_ctg_list, ctg_group_link_dict, ctg_group_dict,
                full_link_dict, linked_ctg_dict, fa_dict, RE_site_dict,
                group_RE_dict, args.max_ctg_len, args.min_RE_sites,
                args.min_links, args.min_link_density, args.min_density_ratio,
                args.ambiguous_cutoff, args.min_group_len, whitelist, n + 1
            )
            if n > 0 and last_round == ctg_group_dict:
                logger.info('[result::round{}] 结果在第 {} 轮后收敛，终止迭代'.format(n + 1, n))
                break
            last_round = ctg_group_dict.copy()

        # 额外救援轮次
        if not args.no_additional_rescue:
            run_reassignment(
                sorted_ctg_list, ctg_group_link_dict, ctg_group_dict,
                full_link_dict, linked_ctg_dict, fa_dict, RE_site_dict,
                group_RE_dict, args.max_ctg_len, args.min_RE_sites,
                args.min_links, args.min_link_density, args.min_density_ratio,
                args.ambiguous_cutoff, args.min_group_len, whitelist, 0
            )

        # 输出重分配结果
        os.mkdir('reassigned_groups')
        group_ctg_dict, group_hiconf_RE_dict = stat_clusters(ctg_group_dict, fa_dict, grouped_ctgs)
        new_ctg_group_dict, new_group_ctg_dict, new_old_group_dict = clusters_output(
            group_ctg_dict, fa_dict, 'reassigned'
        )

        logger.info('重分配完成，耗时 {}s'.format(time.time() - start_time))

        run_hc = False
        # 凝聚层次聚类（可选）
        if args.nclusters and args.nclusters < len(new_group_ctg_dict):
            os.mkdir('hc_groups')
            hc_dict = agglomerative_hierarchical_clustering(
                full_link_dict, grouped_ctgs, new_ctg_group_dict,
                group_hiconf_RE_dict, new_old_group_dict,
                int(args.nclusters), normalize_by_nlinks=args.normalize_by_nlinks
            )
            hc_group_dict = stat_hc_clusters(new_group_ctg_dict, hc_dict)
            new_ctg_group_dict, new_group_ctg_dict, _ = clusters_output(hc_group_dict, fa_dict, 'hc')
            run_hc = True
        elif args.nclusters == 0:
            logger.info('--nclusters 为 0，跳过凝聚层次聚类')
        elif args.nclusters == len(new_group_ctg_dict):
            logger.info(
                '--nclusters（{}）等于当前聚类数（{}），跳过凝聚层次聚类'.format(
                    args.nclusters, len(new_group_ctg_dict))
            )
        else:
            logger.info(
                '--nclusters（{}）大于当前聚类数（{}），建议尝试更高的膨胀值'.format(
                    args.nclusters, len(new_group_ctg_dict))
            )

        # 按组拆分 CLM 文件
        subdir = 'hc_groups' if run_hc else 'reassigned_groups'
        split_clm_file(args.clm, new_group_ctg_dict, new_ctg_group_dict, subdir)

        logger.info('程序完成，总耗时 {}s'.format(time.time() - start_time))
        os.chdir('..')


def main():
    """程序入口：解析参数并启动主流程。"""
    args = parse_arguments()
    run(args, 'ConHiC_reassign.log')


if __name__ == '__main__':
    main()
