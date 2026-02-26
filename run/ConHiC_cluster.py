import sys
import argparse
import os
import logging
import time
import gc

from copy import deepcopy
import pysam
import pickle
from collections import defaultdict, OrderedDict
import random
from itertools import combinations
from array import array
from portion import closed, empty
import gzip
import scipy.sparse as sp
from math import ceil
from numpy import inf, int32, float32, zeros, quantile, arange, power, allclose, median
from numpy import abs as npabs
from numpy import array as ndarray
from numpy.linalg import matrix_power
from scipy.sparse import coo_matrix, dok_matrix, csc_matrix
from scipy.stats import mode
from scipy.optimize import linear_sum_assignment
from sklearn.preprocessing import normalize
from networkx import Graph, find_cliques, connected_components, shortest_path
from decimal import Decimal
import math

try:
    from sparse_dot_mkl import dot_product_mkl
    INTEL_MKL = True
except:
    INTEL_MKL = False

logging.basicConfig(
    format='%(asctime)s <%(filename)s> [%(funcName)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

pysam.set_verbosity(0)


# ─────────────────────────────────────────────────────────────
# 工具函数：随机分箱与链接过滤
# ─────────────────────────────────────────────────────────────

def remove_selected_bins_from_flank_links(flank_link_dict, selected_bins):
    """
    从侧翼链接字典中剔除所有涉及指定 bin 的条目。

    参数
    ----
    flank_link_dict : dict
        格式为 {(bin1, bin2): 强度值, ...}
    selected_bins : list
        需要排除的 bin 名称列表

    返回
    ----
    dict
        过滤后的侧翼链接字典
    """
    excluded = set(selected_bins)
    return {
        (b1, b2): strength
        for (b1, b2), strength in flank_link_dict.items()
        if b1 not in excluded and b2 not in excluded
    }


def select_random_bins(fa_bin_dict, select_ratio=0.25, random_seed=42):
    """
    对每条 contig 的 bin 列表按指定比例随机采样（结果可复现）。

    参数
    ----
    fa_bin_dict : dict
        {contig: [bin1, bin2, ...]}
    select_ratio : float
        采样比例，默认 0.25
    random_seed : int
        随机种子，保证可复现性

    返回
    ----
    list
        被选中的 bin 名称列表
    """
    random.seed(random_seed)
    chosen = []
    min_bins = math.ceil(1 / select_ratio)
    for ctg, bins in fa_bin_dict.items():
        if len(bins) >= min_bins:
            k = max(1, int(len(bins) * select_ratio))
            chosen.extend(random.sample(bins, k))
    return chosen


def generate_multiple_results(fa_bin_dict, select_ratio=0.25, n_groups=10, base_seed=42):
    """
    用不同随机种子生成多组采样结果，每组结果均可独立复现。

    参数
    ----
    fa_bin_dict : dict
        {contig: [bin1, bin2, ...]}
    select_ratio : float
        每组的采样比例
    n_groups : int
        生成的组数
    base_seed : int
        基础随机种子，各组种子依次为 base_seed, base_seed+1, ...

    返回
    ----
    list[list]
        每个元素是一组被选中的 bin 名称
    """
    return [
        select_random_bins(fa_bin_dict, select_ratio, base_seed + i)
        for i in range(n_groups)
    ]


# ─────────────────────────────────────────────────────────────
# 限制酶识别位点相关工具
# ─────────────────────────────────────────────────────────────

def parse_RE_sites(sites):
    """
    将含有模糊碱基 N 的识别位点展开为所有确定性序列。
    该函数会递归展开，直至所有位点均不含 N 为止。

    参数
    ----
    sites : list[str]
        可能含 N 的识别位点列表

    返回
    ----
    list[str]
        不含 N 的所有等价识别位点
    """
    expanded = []
    for site in sites:
        if 'N' in site:
            for base in ('A', 'T', 'C', 'G'):
                expanded.append(site.replace('N', base, 1))
        else:
            expanded.append(site)

    if 'N' not in ''.join(expanded):
        return expanded
    return parse_RE_sites(expanded)


def count_RE_sites(seq, RE):
    """
    统计给定序列中所有限制酶识别位点的总出现次数。

    参数
    ----
    seq : str
        DNA 序列（大写）
    RE : str
        逗号分隔的识别位点字符串

    返回
    ----
    int
        识别位点总计数
    """
    raw_sites = [s.strip().upper() for s in RE.split(',') if s.strip()]
    all_sites = parse_RE_sites(raw_sites)
    return sum(seq.count(s) for s in all_sites)


# ─────────────────────────────────────────────────────────────
# 输入文件解析
# ─────────────────────────────────────────────────────────────

def parse_fasta(fasta, RE='GATC', keep_letter_case=False, logger=logger):
    """
    读取 FASTA 文件并构建 contig 信息字典。

    每条 contig 的值为列表 [序列, 长度, RE位点数+1]，
    其中 RE 位点数额外加 1 以避免后续计算中出现除零错误。

    参数
    ----
    fasta : str
        输入的 FASTA 文件路径
    RE : str
        限制酶识别位点，多个位点用逗号分隔
    keep_letter_case : bool
        是否保留原始字母大小写；默认全部转大写
    logger : logging.Logger
        日志对象

    返回
    ----
    dict
        {contig_name: [sequence, length, RE_count]}
    """
    logger.info('读取输入 FASTA 文件…')

    fa_dict = {}
    with open(fasta) as fh:
        for line in fh:
            if not line.strip():
                continue
            if line.startswith('>'):
                ctg = line.split()[0][1:]
                fa_dict[ctg] = []
            else:
                fa_dict[ctg].append(
                    line.strip() if keep_letter_case else line.strip().upper()
                )

    for ctg, seq_parts in fa_dict.items():
        seq = ''.join(seq_parts)
        re_cnt = count_RE_sites(seq, RE) + 1
        fa_dict[ctg] = [seq, len(seq), re_cnt]

    return fa_dict


def determine_int_type(fa_dict, logger=logger):
    """
    根据 contig 长度分布，判断坐标与距离所需的整型精度。

    当最长 contig 超过 INT32 上限时使用 int64，否则使用 int32；
    距离精度则取决于最长与次长 contig 之和是否超过 INT32 上限。

    参数
    ----
    fa_dict : dict
        parse_fasta 返回的 contig 信息字典
    logger : logging.Logger
        日志对象

    返回
    ----
    tuple[str, str]
        (坐标整型类型, 距离整型类型)，值为 'int32' 或 'int64'
    """
    MAX_INT32 = 2 ** 31 - 1

    lens = sorted([info[1] for info in fa_dict.values()])
    max_len = lens[-1]
    sec_len = 0 if len(fa_dict) < 2 else lens[-2]

    pos_type = 'int64' if max_len > MAX_INT32 else 'int32'
    dist_type = 'int64' if max_len + sec_len > MAX_INT32 else 'int32'

    logger.info(
        '最长 contig：{} bp，次长 contig：{} bp；'
        '坐标精度：{}，距离精度：{}'.format(max_len, sec_len, pos_type, dist_type)
    )

    if pos_type == 'int64':
        logger.warning(
            '存在长度超过 {} bp 的 contig，在 Juicebox 中可视化时可能出现问题'.format(MAX_INT32)
        )

    return pos_type, dist_type


def parse_gfa(gfa_list, fa_dict, logger=logger):
    """
    从一个或多个 GFA 文件中提取 contig 的测序深度和单倍型归属信息。

    参数
    ----
    gfa_list : list[str]
        GFA 文件路径列表
    fa_dict : dict
        parse_fasta 返回的 contig 信息字典（用于长度校验）
    logger : logging.Logger
        日志对象

    返回
    ----
    dict
        {contig_name: (haplotype_index, read_depth)}
    """
    logger.info('读取 GFA 文件…')

    depth_dict = {}
    for idx, gfa in enumerate(gfa_list):
        with open(gfa) as fh:
            for line in fh:
                if not line.startswith('S\t'):
                    continue
                cols = line.split('\t')
                ctg = cols[1]
                gfa_len = int(cols[3].split(':')[-1])
                depth = int(cols[4].split(':')[-1])

                if ctg in fa_dict and gfa_len != fa_dict[ctg][1]:
                    msg = (
                        'contig {} 在 GFA 文件 {} 中的长度与 FASTA 文件不一致，'
                        '请检查输入文件是否匹配'.format(ctg, gfa)
                    )
                    logger.error(msg)
                    raise RuntimeError(msg)

                depth_dict[ctg] = (idx, depth)

    for ctg in fa_dict:
        if ctg not in depth_dict:
            msg = '在 GFA 文件中找不到 contig {}，请检查输入文件'.format(ctg)
            logger.error(msg)
            raise RuntimeError(msg)

    n_gfa, n_fa = len(depth_dict), len(fa_dict)
    if n_gfa > n_fa:
        logger.warning(
            'GFA 文件中的 contig 数量（{}）多于 FASTA 文件（{}），'
            '可能有 contig 在 FASTA 中被移除'.format(n_gfa, n_fa)
        )

    return depth_dict


# ─────────────────────────────────────────────────────────────
# 片段统计与分箱
# ─────────────────────────────────────────────────────────────

def stat_fragments(
        fa_dict, RE, read_depth_dict, whitelist,
        nchrs=0, flank=0, Nx=100, bin_size=0, logger=logger):
    """
    对所有片段（contig 或 bin）进行基础统计，并按需将长 contig 切割为固定大小的 bin。

    当 bin_size > 0 且 contig 长度超过该值时，contig 将被切割为若干 bin；
    否则 contig 作为整体片段保留。

    参数
    ----
    fa_dict : dict
        parse_fasta 返回的 contig 信息字典
    RE : str
        限制酶识别位点
    read_depth_dict : dict
        contig 测序深度字典（可为空）
    whitelist : set
        白名单 contig 集合，其中的片段在过滤步骤中将被强制保留
    nchrs : int
        预期染色体数目，用于自动推算 bin_size
    flank : int
        侧翼区域长度（单位：kbp）；仅统计 contig 两端该范围内的 RE 位点
    Nx : int
        用于筛选代表性片段的 Nx 阈值（0‑100）
    bin_size : int
        分箱大小（单位：bp）；0 表示不分箱，-1 表示自动计算
    logger : logging.Logger
        日志对象

    返回
    ----
    tuple
        (sorted_frag_list, bin_set, bin_size, frag_len_dict,
         Nx_frag_set, RE_site_dict, split_ctg_set, fa_bin_dict, fa_bin_num)
    """

    def count_flank_RE_sites(sequence, length):
        """仅统计片段两端侧翼区域内的 RE 位点数"""
        if not flank or length <= 2 * flank:
            return count_RE_sites(sequence, RE) + 1
        return (
            count_RE_sites(sequence[:flank], RE)
            + count_RE_sites(sequence[length - flank:], RE)
            + 1
        )

    logger.info('统计片段（contig / bin）基础信息…')

    # kbp 转 bp
    flank *= 1000

    total_len = sum(info[1] for info in fa_dict.values())

    if not bin_size:
        logger.info('bin_size 为 {}，不进行分箱'.format(bin_size))
        bin_size = inf
    elif bin_size < 0:
        bin_size = max(min(int(total_len / nchrs / 30), 2_000_000), 100_000)
        logger.info('自动计算 bin_size = {} bp'.format(bin_size))
    else:
        logger.info('手动指定 bin_size = {} bp'.format(bin_size))

    print("bin_size:", bin_size)

    fa_bin_dict, fa_bin_num = {}, {}
    frags = []
    bin_set = set()
    split_ctg_set = set()
    RE_site_dict, frag_len_dict = {}, {}

    for ctg, (seq, ctg_len, RE_sites) in fa_dict.items():
        if ctg_len > bin_size:
            # 长 contig 切割为 bin
            split_ctg_set.add(ctg)
            nbins = ceil(ctg_len / bin_size)
            fa_bin_dict[ctg] = []
            fa_bin_num[ctg] = nbins

            for m in range(nbins):
                bin_name = '{}_bin{}'.format(ctg, m + 1)
                fa_bin_dict[ctg].append(bin_name)
                assert bin_name not in fa_dict
                frags.append(bin_name)
                bin_set.add(bin_name)

                if m + 1 < nbins:
                    blen = bin_size
                    bseq = seq[m * bin_size:(m + 1) * bin_size]
                else:
                    blen = ctg_len - m * bin_size
                    bseq = seq[m * bin_size:]

                RE_site_dict[bin_name] = count_flank_RE_sites(bseq, blen)
                frag_len_dict[bin_name] = blen

                if read_depth_dict:
                    read_depth_dict[bin_name] = read_depth_dict[ctg]

            if read_depth_dict:
                del read_depth_dict[ctg]
        else:
            # 短 contig 整体作为片段
            fa_bin_dict[ctg] = [ctg]
            fa_bin_num[ctg] = 1
            frags.append(ctg)
            frag_len_dict[ctg] = ctg_len

            if not flank or ctg_len <= 2 * flank:
                RE_site_dict[ctg] = RE_sites
            else:
                RE_site_dict[ctg] = count_flank_RE_sites(seq, ctg_len)

        # 序列已不再需要，释放内存
        fa_dict[ctg][0] = None

    # 为避免相同长度片段的排序偏差，先打乱再排序
    random.seed(12345)
    random.shuffle(frags)
    sorted_frag_list = sorted(
        [(f, frag_len_dict[f]) for f in frags],
        key=lambda x: x[1], reverse=True
    )

    # 按 Nx 策略选取代表性片段集合
    len_sum = 0
    Nx_frag_set = set()
    for frag, flen in sorted_frag_list:
        len_sum += flen
        if len_sum / total_len * 100 < Nx or Nx == 100:
            Nx_frag_set.add(frag)

    if Nx != 100:
        Nx_frag_set.add(sorted_frag_list[len(Nx_frag_set)][0])

    # 白名单片段强制加入
    if whitelist:
        for frag, _ in sorted_frag_list:
            if frag.rsplit('_bin', 1)[0] in whitelist:
                Nx_frag_set.add(frag)

    return (
        sorted_frag_list, bin_set, bin_size,
        frag_len_dict, Nx_frag_set, RE_site_dict,
        split_ctg_set, fa_bin_dict, fa_bin_num
    )


# ─────────────────────────────────────────────────────────────
# 辅助判断与矩阵构建
# ─────────────────────────────────────────────────────────────

def is_flank(coord, length, flank):
    """
    判断给定坐标是否位于片段的侧翼区域内。

    当 flank 为 0 时，视所有坐标均在侧翼区域（即使用全长）。

    参数
    ----
    coord : int
        Hi-C 比对坐标（1-based）
    length : int
        片段长度
    flank : int
        侧翼区域长度（bp）

    返回
    ----
    bool
    """
    if not flank:
        return True
    return coord <= flank or coord > length - flank


def dict_to_matrix(link_dict, frag_set, dense_matrix=True, add_self_loops=False):
    """
    将 Hi-C 链接字典转换为对称邻接矩阵。

    参数
    ----
    link_dict : dict
        {(frag_i, frag_j): 链接数}
    frag_set : set
        参与矩阵构建的片段集合
    dense_matrix : bool
        True 返回 ndarray，False 返回 CSC 稀疏矩阵
    add_self_loops : bool
        是否在对角线上添加自环（值为 1），Markov 聚类需要此选项

    返回
    ----
    tuple
        (matrix, frag_index_dict)
        matrix 为邻接矩阵，frag_index_dict 为 {片段名: 行/列索引}
    """
    row, col, data = [], [], []
    index = 0
    shape = len(frag_set)
    frag_index_dict = {}
    seen_frags = set()

    for (fi, fj), links in link_dict.items():
        if fi not in frag_set or fj not in frag_set:
            continue

        seen_frags.add(fi)
        seen_frags.add(fj)

        if fi in frag_index_dict:
            i = frag_index_dict[fi]
        else:
            i = index
            frag_index_dict[fi] = i
            index += 1

        if fj in frag_index_dict:
            j = frag_index_dict[fj]
        else:
            j = index
            frag_index_dict[fj] = j
            index += 1

        row.extend([i, j])
        col.extend([j, i])
        data.extend([links, links])

    # 为没有链接但仍在集合中的片段分配索引
    assert len(seen_frags) == index
    for frag in frag_set - seen_frags:
        frag_index_dict[frag] = index
        index += 1

    if add_self_loops:
        for n in range(shape):
            row.append(n)
            col.append(n)
            data.append(1)

    if not dense_matrix:
        matrix = coo_matrix(
            (data, (row, col)), shape=(shape, shape), dtype=float32
        ).tocsc()
    else:
        matrix = coo_matrix(
            (data, (row, col)), shape=(shape, shape), dtype=float32
        ).toarray()

    return matrix, frag_index_dict


# ─────────────────────────────────────────────────────────────
# CLM 文件与 HT 链接输出
# ─────────────────────────────────────────────────────────────

def output_clm(clm_dict):
    """将 CLM 字典写入 paired_links.clm 文件，供后续排序步骤使用。"""

    logger.info('将 clm_dict 写入 paired_links.clm…')

    ori_tuple = (('+', '+'), ('+', '-'), ('-', '+'), ('-', '-'))

    with open('paired_links.clm', 'w') as fout:
        for (c1, c2), lst in clm_dict.items():
            if len(lst) < 8:
                continue
            for n in range(4):
                vals = ['{0} {0}'.format(v) for v in sorted(lst[n::4])]
                fout.write('{}{} {}{}\t{}\t{}\n'.format(
                    c1, ori_tuple[n][0], c2, ori_tuple[n][1],
                    len(vals) * 2, ' '.join(vals)
                ))


def update_clm_dict(clm_dict, ctg_name_pair, len_i, len_j, coord_i_0, coord_j_0):
    """向 CLM 字典中追加一对 Hi-C 链接的四种方向距离。"""
    clm_dict[ctg_name_pair].extend((
        len_i - coord_i_0 + coord_j_0,
        len_i - coord_i_0 + len_j - coord_j_0,
        coord_i_0 + coord_j_0,
        coord_i_0 + len_j - coord_j_0,
    ))


def update_HT_link_dict(HT_link_dict, ctg_i, ctg_j, len_i, len_j, coord_i, coord_j):
    """
    根据 Hi-C 链接坐标判断比对端点位于 contig 的头（H）还是尾（T），
    并更新半 contig 链接字典。
    """

    def add_suffix(ctg, ctg_len, coord):
        return ctg + ('_T' if coord * 2 > ctg_len else '_H')

    ht_i = add_suffix(ctg_i, len_i, coord_i)
    ht_j = add_suffix(ctg_j, len_j, coord_j)
    HT_link_dict[(ht_i, ht_j)] += 1


# ─────────────────────────────────────────────────────────────
# 等位链接检测与过滤
# ─────────────────────────────────────────────────────────────

def cal_concordance_ratio(coord_list, shorter_len, nwindows):
    """
    通过计算 Hi-C 链接坐标的一致性比率来判断两条 contig 是否为等位序列。

    等位 contig 的 Hi-C 链接通常呈现 y = x + b 或 y = -x + b 的线性模式，
    此函数通过众数频率来量化这种一致性。

    参数
    ----
    coord_list : array
        交错存储的坐标对 [x0, y0, x1, y1, ...]
    shorter_len : int
        两条 contig 中较短者的长度
    nwindows : int
        用于定义窗口宽度的分箱数

    返回
    ----
    float
        0~1 之间的一致性比率，越高越可能为等位关系
    """
    bw = shorter_len // nwindows
    npairs = len(coord_list) // 2
    yx_diff = [(coord_list[2*n+1] - coord_list[2*n]) // bw for n in range(npairs)]
    yx_sum  = [(coord_list[2*n+1] + coord_list[2*n]) // bw for n in range(npairs)]
    return max(
        mode(yx_diff, keepdims=False)[1] / npairs,
        mode(yx_sum,  keepdims=False)[1] / npairs,
    )


def cal_concentration_adj_ratio(coord_list, bin_width=10000):
    """
    计算链接集中度校正系数，用于识别并惩罚高度集中于局部区域的 Hi-C 信号。

    参数
    ----
    coord_list : array
        交错存储的坐标对
    bin_width : int
        统计链接集中度的箱宽（bp）

    返回
    ----
    float
        0~1 之间的校正系数，越接近 1 表示链接分布越均匀
    """
    npairs = len(coord_list) // 2
    x_bins, y_bins = defaultdict(int), defaultdict(int)

    for n in range(npairs):
        x_bins[coord_list[2*n] // bin_width] += 1
        y_bins[coord_list[2*n+1] // bin_width] += 1

    xm = median([v for v in x_bins.values() if v])
    ym = median([v for v in y_bins.values() if v])

    cx = sum(v for v in x_bins.values() if v >= 10 * xm) / npairs
    cy = sum(v for v in y_bins.values() if v >= 10 * ym) / npairs

    return (1 - cx) * (1 - cy)


def record_coord_pairs(ctg_coord_dict, ctg_name_pair, coord_i, coord_j, max_read_pairs, fa_dict, args):
    """
    缓存 Hi-C 坐标对，并在积累足够数量后计算一致性比率或集中度系数。

    参数
    ----
    ctg_coord_dict : dict
        坐标缓存字典；当坐标对达到阈值时，将被替换为计算结果列表
    ctg_name_pair : tuple
        contig 名称对
    coord_i, coord_j : int
        当前 Hi-C 链接的两端坐标
    max_read_pairs : int
        触发计算的最大坐标对数量
    fa_dict : dict
        contig 信息字典
    args : argparse.Namespace
        运行参数
    """
    if not isinstance(ctg_coord_dict[ctg_name_pair], list):
        ctg_coord_dict[ctg_name_pair].extend((coord_i, coord_j))

        if len(ctg_coord_dict[ctg_name_pair]) >= max_read_pairs * 2:
            if args.remove_allelic_links:
                shorter = min(fa_dict[ctg_name_pair[0]][1], fa_dict[ctg_name_pair[1]][1])
                cr = cal_concordance_ratio(ctg_coord_dict[ctg_name_pair], shorter, args.nwindows)
                ctg_coord_dict[ctg_name_pair] = [cr, 1]
            if args.remove_concentrated_links:
                adj = cal_concentration_adj_ratio(ctg_coord_dict[ctg_name_pair])
                if args.remove_allelic_links:
                    ctg_coord_dict[ctg_name_pair][1] = adj
                else:
                    ctg_coord_dict[ctg_name_pair] = [0, adj]


def remove_allelic_HiC_links(
        fa_dict, ctg_coord_dict, full_link_dict, args,
        flank_link_dict=None, filtered_frags=None, ctg_pair_to_frag=None,
        logger=logger):
    """
    识别并移除等位 contig 之间的 Hi-C 链接。

    本函数分两阶段工作：
    （1）基于一致性比率阈值，直接剔除已识别的等位 contig 对之间的链接；
    （2）构建等位 contig 组，通过匈牙利算法（最大二分匹配）确定最优匹配，
         并移除非最优匹配 contig 对之间的链接。

    参数
    ----
    fa_dict : dict
        contig 信息字典
    ctg_coord_dict : dict
        坐标缓存或计算结果字典
    full_link_dict : dict
        完整的 contig 间链接字典（将被原地修改）
    args : argparse.Namespace
        运行参数（包含 ploidy、阈值等）
    flank_link_dict : dict, optional
        侧翼链接字典（将被原地修改）
    filtered_frags : set, optional
        当前过滤后的片段集合
    ctg_pair_to_frag : dict, optional
        contig 对到 bin 对的映射（分箱模式下使用）
    logger : logging.Logger
        日志对象

    返回
    ----
    set or None
        过滤等位链接后剩余的片段集合（仅在提供 flank_link_dict 时返回）
    """

    def update_link_dicts(ctg_name_pair, link_type):
        if link_type == 1:
            inter_allele_dict[ctg_name_pair] = full_link_dict[ctg_name_pair]
            allelic_ctg_set.add(ctg_name_pair[0])
            allelic_ctg_set.add(ctg_name_pair[1])

        del full_link_dict[ctg_name_pair]

        if flank_link_dict:
            if ctg_pair_to_frag:
                for fp in ctg_pair_to_frag[ctg_name_pair]:
                    if fp in flank_link_dict and fp[0] in filtered_frags and fp[1] in filtered_frags:
                        del flank_link_dict[fp]
            elif (ctg_name_pair in flank_link_dict
                  and ctg_name_pair[0] in filtered_frags
                  and ctg_name_pair[1] in filtered_frags):
                del flank_link_dict[ctg_name_pair]

    def get_weakest_edge(graph):
        weakest = (None, None, inf)
        for n1, n2, d in graph.edges(data=True):
            if n1 != n2 and d['weight'] < weakest[-1]:
                weakest = (n1, n2, d['weight'])
        assert weakest[0] is not None
        return weakest

    def split_cliques(graph, cliques, ploidy, cached):
        new_cliques = set()
        for clique in cliques:
            clique = tuple(clique)
            if len(clique) > ploidy:
                if clique not in cached:
                    sub = Graph(graph.subgraph(clique))
                    n1, n2, _ = get_weakest_edge(sub)
                    sub.remove_edge(n1, n2)
                    cached.add(clique)
                    new_cliques |= split_cliques(sub, find_cliques(sub), ploidy, cached)
            else:
                new_cliques.add(tuple(clique))
        return new_cliques

    def solve_max_matching(group_pair):
        g1, g2 = group_pair
        degree = max(len(g1), len(g2))
        mat = zeros((degree, degree), dtype=int)
        for i1, c1 in enumerate(g1):
            for i2, c2 in enumerate(g2):
                pair = tuple(sorted([c1, c2]))
                if pair in full_link_dict:
                    mat[i1, i2] = full_link_dict[pair]
        return linear_sum_assignment(-mat)

    ploidy = args.remove_allelic_links
    min_pairs = args.min_read_pairs
    cr_cutoff = args.concordance_ratio_cutoff

    inter_allele_dict = {}
    allelic_ctg_set = set()

    # 阶段一：基于一致性比率剔除等位链接
    for pair, data in ctg_coord_dict.items():
        if isinstance(data, list):
            logger.debug('{} {} links={} cr={}'.format(*pair, full_link_dict[pair], data[0]))
            if data[0] > cr_cutoff:
                update_link_dicts(pair, 1)
        else:
            if len(data) >= min_pairs * 2:
                shorter = min(fa_dict[pair[0]][1], fa_dict[pair[1]][1])
                cr = cal_concordance_ratio(data, shorter, args.nwindows)
                logger.debug('{} {} links={} cr={}'.format(*pair, full_link_dict[pair], cr))
                if cr > cr_cutoff:
                    update_link_dicts(pair, 1)
            else:
                logger.debug('{} {} links={} cr=0'.format(*pair, full_link_dict[pair]))

    # 阶段二（多倍体）：通过 clique 分割和最大匹配进一步过滤
    if ploidy > 2:
        allele_matrix, ctg_idx = dict_to_matrix(inter_allele_dict, allelic_ctg_set)
        idx_ctg = {i: c for c, i in ctg_idx.items()}
        g = Graph(allele_matrix)
        allele_groups = split_cliques(g, find_cliques(g), ploidy, set())
        gc.collect()
        unique_groups = set()
        for grp in allele_groups:
            unique_groups.add(tuple(sorted([idx_ctg[i] for i in grp])))
    else:
        unique_groups = set(inter_allele_dict.keys())

    ctg_allele_group = defaultdict(set)
    for grp in unique_groups:
        for c in grp:
            ctg_allele_group[c].add(grp)

    solution_cache = {}
    nonmax_pairs = set()

    for pair in full_link_dict:
        c1, c2 = pair
        if c1 not in ctg_allele_group or c2 not in ctg_allele_group:
            continue

        for g1 in ctg_allele_group[c1]:
            for g2 in ctg_allele_group[c2]:
                gp = tuple(sorted([g1, g2]))
                if gp not in solution_cache:
                    solution_cache[gp] = solve_max_matching(gp)
                sol = solution_cache[gp]

                if c1 in gp[0]:
                    i1 = gp[0].index(c1)
                    i2 = gp[1].index(c2)
                else:
                    i1 = gp[0].index(c2)
                    i2 = gp[1].index(c1)

                if sol[1][i1] != i2:
                    nonmax_pairs.add(pair)
                    break
            else:
                continue
            break

    for pair in nonmax_pairs:
        logger.debug('{} {} links={} non-maximum matching'.format(*pair, full_link_dict[pair]))
        update_link_dicts(pair, 2)

    if flank_link_dict:
        remaining = set()
        for f1, f2 in flank_link_dict:
            if f1 in filtered_frags and f2 in filtered_frags:
                remaining.add(f1)
                remaining.add(f2)

        removed = filtered_frags - remaining
        logger.info('过滤等位链接后移除孤立片段：移除 {}，保留 {}'.format(len(removed), len(remaining)))
        for frag in removed:
            logger.debug('片段 {} 孤立，已移除'.format(frag))

        return remaining


# ─────────────────────────────────────────────────────────────
# 单倍型间链接降权
# ─────────────────────────────────────────────────────────────

def reduce_inter_hap_HiC_links(link_dict, read_depth_dict, phasing_weight, target='flank_link_dict'):
    """
    根据单倍型归属信息，对跨单倍型的 Hi-C 链接进行降权处理。

    参数
    ----
    link_dict : dict
        待修改的链接字典（原地修改）
    read_depth_dict : dict
        {片段名: (单倍型索引, 测序深度)}
    phasing_weight : float
        降权比例；1.0 表示完全移除跨单倍型链接
    target : str
        日志标识字符串
    """
    logger.info('对 {} 中的跨单倍型链接进行降权…'.format(target))

    to_remove = set()
    for pair in link_dict:
        if read_depth_dict[pair[0]][0] != read_depth_dict[pair[1]][0]:
            link_dict[pair] -= link_dict[pair] * phasing_weight
            if link_dict[pair] == 0:
                to_remove.add(pair)

    for pair in to_remove:
        del link_dict[pair]


# ─────────────────────────────────────────────────────────────
# 序列化工具
# ─────────────────────────────────────────────────────────────

def output_pickle(obj, from_, to):
    """将任意 Python 对象序列化为 pickle 文件。"""
    logger.info('将 {} 写入 {}…'.format(from_, to))
    with open(to, 'wb') as fh:
        pickle.dump(obj, fh)


def load_pickle(path):
    """从 pickle 文件中反序列化并返回对象。"""
    with open(path, 'rb') as fh:
        return pickle.load(fh)


# ─────────────────────────────────────────────────────────────
# 链接归一化
# ─────────────────────────────────────────────────────────────

def normalize_by_nlinks(flank_link_dict, frag_link_dict):
    """
    以两端片段各自总链接数的几何平均值对侧翼链接进行归一化。

    参数
    ----
    flank_link_dict : dict
        侧翼链接字典（原地修改）
    frag_link_dict : dict
        {片段名: 总链接数}
    """
    logger.info('按总链接数对侧翼链接进行归一化…')
    for fi, fj in flank_link_dict:
        flank_link_dict[(fi, fj)] /= (frag_link_dict[fi] * frag_link_dict[fj]) ** 0.5


def normalize_by_length(flank_link_dict, frag_len_dict, flank):
    """
    以两端片段侧翼区域长度之积对链接数进行归一化。

    参数
    ----
    flank_link_dict : dict
        侧翼链接字典（原地修改）
    frag_len_dict : dict
        {片段名: 片段长度（bp）}
    flank : int
        侧翼区域长度（kbp）
    """
    logger.info('按片段长度对侧翼链接进行归一化…')
    two_flanks = flank * 2000
    for fi, fj in flank_link_dict:
        li = min(frag_len_dict[fi], two_flanks)
        lj = min(frag_len_dict[fj], two_flanks)
        flank_link_dict[(fi, fj)] /= (li / 1_000_000) * (lj / 1_000_000)


# ─────────────────────────────────────────────────────────────
# 片段过滤
# ─────────────────────────────────────────────────────────────

def filter_fragments(
        Nx_frag_set, RE_site_dict, RE_site_cutoff, frag_link_dict, density_lower,
        density_upper, topN, rank_sum_upper, rank_sum_hard_cutoff, flank_link_dict,
        read_depth_dict, read_depth_upper, whitelist):
    """
    对片段依次进行多轮过滤，剔除低质量、高噪声和嵌合片段。

    过滤流程：
    （1）Nx 过滤：仅保留长度达到 Nx 的片段；
    （2）RE 位点过滤：移除 RE 位点数低于阈值的片段；
    （3）链接密度过滤：移除密度极低或极高的片段；
    （4）测序深度过滤（可选）：移除深度异常偏高的片段；
    （5）秩和过滤：移除邻域链接分布异常的片段。

    参数
    ----
    Nx_frag_set : set
        通过 Nx 筛选的片段集合
    RE_site_dict : dict
        {片段名: RE 位点数}
    RE_site_cutoff : int
        RE 位点数下限
    frag_link_dict : dict
        {片段名: 总链接数}
    density_lower, density_upper : str
        链接密度下/上限，支持分数模式和倍数模式（以 X 结尾）
    topN : int
        秩和计算时使用的最近邻片段数
    rank_sum_upper : str
        秩和上限
    rank_sum_hard_cutoff : int
        秩和硬过滤阈值（0 表示禁用）
    flank_link_dict : dict
        侧翼链接字典，用于构建邻接矩阵
    read_depth_dict : dict
        测序深度字典（可为空以跳过深度过滤）
    read_depth_upper : str
        深度上限
    whitelist : set
        白名单 contig 集合

    返回
    ----
    set
        过滤后保留的片段集合
    """
    logger.info('过滤片段…')

    frags_in_whitelist = set()
    density_list = []
    total_links = 0
    total_RE = 1

    # 步骤 1 & 2：Nx + RE 位点过滤
    for frag in Nx_frag_set:
        re_cnt = RE_site_dict[frag]
        if re_cnt > RE_site_cutoff:
            if frag in frag_link_dict:
                links = frag_link_dict[frag]
                total_links += links
                total_RE += re_cnt - 1
                density_list.append((frag, links / re_cnt))
            else:
                density_list.append((frag, 0))
        if whitelist and frag.rsplit('_bin', 1)[0] in whitelist:
            frags_in_whitelist.add(frag)

    nx_n = len(Nx_frag_set)
    logger.info('[Nx 过滤] 保留 {} 个片段'.format(nx_n))
    logger.info('[RE 位点过滤] 移除 {}，保留 {}'.format(nx_n - len(density_list), len(density_list)))

    # 步骤 3：链接密度过滤
    density_list.sort(key=lambda x: x[1])
    avg_density = total_links / total_RE
    n_remain = len(density_list)

    p_lower = check_param('--density_lower', density_lower, {'X', 'x'})
    p_upper = check_param('--density_upper', density_upper, {'X', 'x'})

    if p_lower[-1] in {'X', 'x'}:
        for lower, (_, d) in enumerate(density_list):
            if d >= avg_density * p_lower[0]:
                break
        else:
            lower += 1
        logger.info('[密度过滤] --density_lower {} 等效分数模式：{}'.format(
            density_lower, lower / n_remain))
    else:
        lower = int(n_remain * float(density_lower))

    if p_upper[-1] in {'X', 'x'}:
        for upper, (_, d) in enumerate(density_list):
            if d > avg_density * p_upper[0]:
                break
        else:
            upper += 1
        logger.info('[密度过滤] --density_upper {} 等效分数模式：{}'.format(
            density_upper, upper / n_remain))
    else:
        upper = int(n_remain * float(density_upper))

    filtered_frags = {f for f, _ in density_list[lower:upper]}
    logger.info('[密度过滤] 移除 {}，保留 {}'.format(
        n_remain - len(filtered_frags), len(filtered_frags)))

    density_list_full = density_list
    density_list = density_list[lower:upper]

    # 步骤 4：测序深度过滤（可选）
    if read_depth_dict:
        depth_list = sorted(
            [(f, read_depth_dict[f][1]) for f, _ in density_list_full],
            key=lambda x: x[1]
        )
        p_depth = check_param('--read_depth_upper', read_depth_upper, {'X', 'x'})
        q1, m, q3 = quantile([d for _, d in depth_list], (0.25, 0.5, 0.75))
        iqr = q3 - q1
        logger.info('[深度过滤] Q1={}, 中位数={}, Q3={}, IQR={}'.format(q1, m, q3, iqr))

        if p_depth[-1]:
            limit = q3 + p_depth[0] * iqr
            for upper, (_, d) in enumerate(depth_list):
                if d > limit:
                    break
            else:
                upper += 1
        else:
            upper = int(n_remain * float(read_depth_upper))

        filtered_frags &= {f for f, _ in depth_list[:upper]}
        removed_by_depth = {f for f, _ in depth_list[upper:]}
        removed_by_density = {f for f, _ in density_list_full[:lower] + density_list_full[upper:]}
        specific = removed_by_depth - removed_by_density
        logger.info('[深度过滤] 移除 {}，保留 {}'.format(len(specific), len(filtered_frags)))

        density_list = [(f, d) for f, d in density_list if f in filtered_frags]

    # 步骤 5：秩和过滤
    mat, frag_idx = dict_to_matrix(flank_link_dict, filtered_frags)
    idx_frag = {i: f for f, i in frag_idx.items()}

    rank_dict = {}
    for frag, _ in density_list:
        idx = frag_idx[frag]
        ranked = sorted(enumerate(mat[idx, :]), key=lambda x: x[1], reverse=True)
        rank_dict[frag] = [idx_frag[i] for i, _ in ranked]

    rank_sum_list = []
    hard_removed = 0
    for frag, _ in density_list:
        topn = rank_dict[frag][:topN]
        rs = sum(
            min(rank_dict[f1].index(f2), rank_dict[f2].index(f1))
            for f1, f2 in combinations(topn, 2)
        )
        if rank_sum_hard_cutoff and rs > rank_sum_hard_cutoff:
            hard_removed += 1
            logger.debug('[秩和过滤] 片段 {} 被硬过滤，秩和={}'.format(frag, rs))
            continue
        rank_sum_list.append((frag, rs))

    rank_sum_list.sort(key=lambda x: x[1])
    n_rs = len(rank_sum_list)

    if rank_sum_hard_cutoff:
        logger.info('[秩和过滤] 硬过滤移除 {}，保留 {}'.format(hard_removed, n_rs))

    p_rs = check_param('--rank_sum_upper', rank_sum_upper, {'X', 'x'})
    q1, m, q3 = quantile([rs for _, rs in rank_sum_list], (0.25, 0.5, 0.75))
    iqr = q3 - q1
    logger.info('[秩和过滤] Q1={}, 中位数={}, Q3={}, IQR={}'.format(q1, m, q3, iqr))

    if p_rs[-1]:
        limit = q3 + p_rs[0] * iqr
        for upper, (_, rs) in enumerate(rank_sum_list):
            if rs > limit:
                break
        else:
            upper += 1
    else:
        upper = int(n_rs * float(rank_sum_upper))

    filtered_frags = {f for f, _ in rank_sum_list[:upper]}
    logger.info('[秩和过滤] 移除 {}，保留 {}'.format(n_rs - len(filtered_frags), len(filtered_frags)))

    # 补充白名单片段
    if frags_in_whitelist:
        added = 0
        for frag in frags_in_whitelist:
            if frag not in filtered_frags:
                added += 1
                filtered_frags.add(frag)
                logger.debug('[秩和过滤] 白名单片段 {} 被强制加入'.format(frag))
        logger.info('[秩和过滤] 白名单补充 {} 个，最终用于聚类：{} 个'.format(added, len(filtered_frags)))

    return filtered_frags


# ─────────────────────────────────────────────────────────────
# 比对文件解析
# ─────────────────────────────────────────────────────────────

def check_sorting_order(fbam):
    """检查 BAM 文件的排序方式，若为坐标排序则报错退出。"""
    hd = fbam.header.get('HD')
    if hd is not None and 'SO' in hd:
        so = hd['SO']
        if so in {'unsorted', 'queryname'}:
            logger.info('BAM 文件排序方式：{}'.format(so))
        elif so == 'coordinate':
            msg = 'BAM 文件为坐标排序，应为未排序或按 read 名排序'.format(so)
            logger.error(msg)
            raise RuntimeError(msg)
    else:
        logger.warning('无法确定 BAM 文件排序方式，程序将继续运行')


def parse_alignments_for_ctgs(
        alignments, fa_dict, args, ctg_len_dict, Nx_ctg_set, pos_int_type, dist_int_type):
    """
    在不进行分箱的情况下，解析 Hi-C 比对并构建各类链接字典。

    参数
    ----
    alignments : iterable
        Hi-C 比对记录，每条为 (ref, mref, pos, mpos)
    fa_dict : dict
        contig 信息字典
    args : argparse.Namespace
        运行参数
    ctg_len_dict : dict
        {contig: 长度}
    Nx_ctg_set : set
        通过 Nx 筛选的 contig 集合
    pos_int_type : str
        坐标整型类型
    dist_int_type : str
        距离整型类型

    返回
    ----
    tuple
        (full_link_dict, flank_link_dict, HT_link_dict,
         clm_dict, ctg_link_dict, ctg_coord_dict)
    """
    logger.info('解析比对文件…')

    flank = args.flank * 1000

    full_link_dict = defaultdict(int)
    ctg_coord_dict = defaultdict(
        lambda: array('i') if pos_int_type == 'int32' else array('l')
    )
    flank_link_dict = defaultdict(int)
    HT_link_dict = defaultdict(int)
    ctg_link_dict = defaultdict(int)
    clm_dict = defaultdict(
        lambda: array('i') if dist_int_type == 'int32' else array('l')
    )

    for ref, mref, pos, mpos in alignments:
        if ref not in fa_dict or mref not in fa_dict:
            continue

        (ci, xi), (cj, xj) = sorted(((ref, pos + 1), (mref, mpos + 1)))
        pair = (ci, cj)
        li, lj = ctg_len_dict[ci], ctg_len_dict[cj]

        if ci in Nx_ctg_set and cj in Nx_ctg_set and is_flank(xi, li, flank) and is_flank(xj, lj, flank):
            flank_link_dict[pair] += 1
            ctg_link_dict[ci] += 1
            ctg_link_dict[cj] += 1

        update_clm_dict(clm_dict, pair, li, lj, xi - 1, xj - 1)
        update_HT_link_dict(HT_link_dict, ci, cj, li, lj, xi, xj)
        full_link_dict[pair] += 1

        if args.remove_allelic_links or args.remove_concentrated_links:
            record_coord_pairs(ctg_coord_dict, pair, xi, xj, args.max_read_pairs, fa_dict, args)

    return full_link_dict, flank_link_dict, HT_link_dict, clm_dict, ctg_link_dict, ctg_coord_dict


def parse_alignments(
        alignments, fa_dict, args, bin_size, frag_len_dict,
        Nx_frag_set, split_ctg_set, pos_int_type, dist_int_type):
    """
    在部分 contig 被切割为 bin 的情况下，解析 Hi-C 比对。

    参数
    ----
    alignments : iterable
        Hi-C 比对记录，每条为 (ref, mref, pos, mpos)
    fa_dict : dict
        contig 信息字典
    args : argparse.Namespace
        运行参数
    bin_size : int / float
        分箱大小（bp 或 inf）
    frag_len_dict : dict
        {片段名: 长度}
    Nx_frag_set : set
        通过 Nx 筛选的片段集合
    split_ctg_set : set
        被切割的 contig 集合
    pos_int_type : str
        坐标整型类型
    dist_int_type : str
        距离整型类型

    返回
    ----
    tuple
        (full_link_dict, flank_link_dict, HT_link_dict, clm_dict,
         frag_link_dict, ctg_coord_dict, ctg_pair_to_frag)
    """

    def convert_frags(ctg, coord):
        if ctg in split_ctg_set:
            nb = ceil(coord / bin_size)
            return '{}_bin{}'.format(ctg, nb), coord - (nb - 1) * bin_size, True
        return ctg, coord, False

    logger.info('解析比对文件（含分箱）…')

    flank = args.flank * 1000
    full_link_dict = defaultdict(int)
    ctg_coord_dict = defaultdict(
        lambda: array('i') if pos_int_type == 'int32' else array('l')
    )
    ctg_pair_to_frag = defaultdict(set)
    flank_link_dict = defaultdict(int)
    HT_link_dict = defaultdict(int)
    frag_link_dict = defaultdict(int)
    clm_dict = defaultdict(
        lambda: array('i') if dist_int_type == 'int32' else array('l')
    )

    for ref, mref, pos, mpos in alignments:
        if ref == mref and ref not in split_ctg_set:
            continue
        if ref not in fa_dict or mref not in fa_dict:
            continue

        (ci, xi), (cj, xj) = sorted(((ref, pos + 1), (mref, mpos + 1)))
        ctg_pair = (ci, cj)

        fi, fxi, i_bin = convert_frags(ci, xi)
        fj, fxj, j_bin = convert_frags(cj, xj)

        if fi == fj:
            continue

        if i_bin or j_bin:
            (fi, fxi), (fj, fxj) = sorted(((fi, fxi), (fj, fxj)))

        fpair = (fi, fj)
        fli, flj = frag_len_dict[fi], frag_len_dict[fj]

        if fi in Nx_frag_set and fj in Nx_frag_set and is_flank(fxi, fli, flank) and is_flank(fxj, flj, flank):
            flank_link_dict[fpair] += 1
            frag_link_dict[fi] += 1
            frag_link_dict[fj] += 1

        if args.remove_allelic_links:
            ctg_pair_to_frag[ctg_pair].add(fpair)

        if ref != mref:
            full_link_dict[ctg_pair] += 1
            li, lj = fa_dict[ci][1], fa_dict[cj][1]
            update_clm_dict(clm_dict, ctg_pair, li, lj, xi - 1, xj - 1)
            update_HT_link_dict(HT_link_dict, ci, cj, li, lj, xi, xj)

            if args.remove_allelic_links or args.remove_concentrated_links:
                record_coord_pairs(ctg_coord_dict, ctg_pair, xi, xj, args.max_read_pairs, fa_dict, args)

    return full_link_dict, flank_link_dict, HT_link_dict, clm_dict, frag_link_dict, ctg_coord_dict, ctg_pair_to_frag


# ─────────────────────────────────────────────────────────────
# 比对记录生成器
# ─────────────────────────────────────────────────────────────

def pairs_generator(pairs, aln_format):
    """逐行读取 pairs 文件并生成 (ref, mref, pos, mpos) 四元组。"""
    fopen = open if aln_format == 'pairs' else gzip.open
    with fopen(pairs, 'rt') as fh, open('alignments.bed', 'w') as bed:
        for line in fh:
            if not line.strip() or line.startswith('#'):
                continue
            cols = line.split()
            ref, pos, mref, mpos = cols[1], int(cols[2]) - 1, cols[3], int(cols[4]) - 1
            bed.write(
                '{0}\t{1}\t{2}\t{3}/1\t255\t.\n'
                '{4}\t{5}\t{6}\t{3}/2\t255\t.\n'.format(
                    ref, pos, pos, cols[0], mref, mpos, mpos)
            )
            yield ref, mref, pos, mpos


def pairs_generator_inter_ctgs(pairs, aln_format):
    """逐行读取 pairs 文件，仅生成跨 contig 的 Hi-C 记录。"""
    fopen = open if aln_format == 'pairs' else gzip.open
    with fopen(pairs, 'rt') as fh, open('alignments.bed', 'w') as bed:
        for line in fh:
            if not line.strip() or line.startswith('#'):
                continue
            cols = line.split()
            ref, pos, mref, mpos = cols[1], int(cols[2]) - 1, cols[3], int(cols[4]) - 1
            bed.write(
                '{0}\t{1}\t{2}\t{3}/1\t255\t.\n'
                '{4}\t{5}\t{6}\t{3}/2\t255\t.\n'.format(
                    ref, pos, pos, cols[0], mref, mpos, mpos)
            )
            if ref != mref:
                yield ref, mref, pos, mpos


def bam_generator(bam, threads, format_options):
    """逐条读取 BAM 文件并生成 (ref, mref, pos, mpos) 四元组。"""
    with pysam.AlignmentFile(bam, mode='rb', threads=threads, format_options=format_options) as fh:
        check_sorting_order(fh)
        for aln in fh:
            yield aln.reference_name, aln.next_reference_name, aln.reference_start, aln.next_reference_start


# ─────────────────────────────────────────────────────────────
# Markov 聚类
# ─────────────────────────────────────────────────────────────

def prune(matrix, pruning, dense_matrix):
    """
    对 Markov 聚类中间矩阵进行剪枝，保留每列最大值并将小值清零。

    参数
    ----
    matrix : array or sparse matrix
        当前迭代的转移概率矩阵
    pruning : float
        小于此阈值的元素将被置零
    dense_matrix : bool
        是否为密集矩阵

    返回
    ----
    array or sparse matrix
        剪枝并重新归一化后的矩阵
    """
    if not dense_matrix and matrix.nnz / matrix.shape[0] ** 2 < 0.05:
        pruned = dok_matrix(matrix.shape, dtype=float32)
        mask = matrix >= pruning
        pruned[mask] = matrix[mask]
        pruned = pruned.tocsc()
    else:
        pruned = matrix.toarray() if not dense_matrix else matrix.copy()
        pruned[pruned < pruning] = 0
        if not dense_matrix:
            pruned = csc_matrix(pruned)

    ncols = matrix.shape[1]
    cols = arange(ncols)
    rows = matrix.argmax(axis=0).reshape((ncols,))
    pruned[rows, cols] = matrix[rows, cols]

    return normalize(pruned, norm='l1', axis=0)


def mkl_matrix_power(matrix, n):
    """使用 Intel MKL 递归计算矩阵的 n 次幂。"""
    if n == 2:
        return dot_product_mkl(matrix, matrix)
    return dot_product_mkl(matrix, mkl_matrix_power(matrix, n - 1))


def mcl(matrix, expansion, inflation, iters, pruning, dense_matrix):
    """
    执行 Markov 聚类算法的核心迭代过程。

    迭代步骤：展开（矩阵乘幂）→ 膨胀（逐元素乘幂）→ 剪枝 → 收敛判断。

    参数
    ----
    matrix : array or sparse matrix
        初始概率矩阵
    expansion : int
        展开参数（矩阵幂次）
    inflation : float
        膨胀参数
    iters : int
        最大迭代次数
    pruning : float
        剪枝阈值
    dense_matrix : bool
        是否使用密集矩阵

    返回
    ----
    array or sparse matrix
        收敛后的矩阵
    """
    last_matrix = None

    for n in range(iters):
        if n != 0:
            if not dense_matrix:
                matrix = mkl_matrix_power(matrix, expansion)
            else:
                matrix = matrix_power(matrix, expansion)

        if not dense_matrix:
            matrix = normalize(matrix.power(inflation), norm='l1', axis=0)
        else:
            matrix = normalize(power(matrix, inflation), norm='l1', axis=0)

        matrix = prune(matrix, pruning, dense_matrix)

        if n > 1:
            if not dense_matrix:
                d = npabs(matrix - last_matrix) - 1e-5 * npabs(last_matrix)
                if d.max() <= 1e-8:
                    logger.info(
                        '矩阵在第 {} 轮迭代后收敛（展开={}，膨胀={}，最大迭代={}，剪枝={}）'.format(
                            n + 1, expansion, inflation, iters, pruning)
                    )
                    return matrix
            elif allclose(matrix, last_matrix):
                logger.info(
                    '矩阵在第 {} 轮迭代后收敛（展开={}，膨胀={}，最大迭代={}，剪枝={}）'.format(
                        n + 1, expansion, inflation, iters, pruning)
                )
                return matrix

        last_matrix = matrix.copy()

    logger.info(
        '矩阵在 {} 轮迭代后未收敛（展开={}，膨胀={}，最大迭代={}，剪枝={}）'.format(
            n + 1, expansion, inflation, iters, pruning)
    )
    return matrix


def interpret_result(result_matrix, dense_matrix):
    """
    从收敛后的 Markov 矩阵中提取聚类结果。

    吸引子（attractor）是对角线上非零的行，每个吸引子对应一个聚类。

    参数
    ----
    result_matrix : array or sparse matrix
        收敛后的矩阵
    dense_matrix : bool
        是否为密集矩阵

    返回
    ----
    list or None
        聚类列表，每个元素为列索引元组；若结果不完整则返回 None
    """
    shape = result_matrix.shape[0]
    attractors = result_matrix.diagonal().nonzero()[0]
    clusters = set()

    for att in attractors:
        if not dense_matrix:
            cluster = tuple(result_matrix.getrow(att).nonzero()[1].tolist())
        else:
            cluster = tuple(n for n in range(shape) if result_matrix[att, n] != 0)
        clusters.add(cluster)

    nodes = set()
    for c in clusters:
        for n in c:
            if n in nodes:
                return None
            nodes.add(n)

    if len(nodes) != shape:
        return None

    return list(clusters)


def get_main_groups(result_clusters, len_ratio):
    """
    基于相邻聚类长度比值确定主要聚类数目。

    参数
    ----
    result_clusters : list
        按长度降序排列的聚类列表
    len_ratio : float
        相邻聚类长度比值阈值

    返回
    ----
    int
        主要聚类数目
    """
    len_sum = 0
    main_groups = len(result_clusters)
    for i in range(len(result_clusters) - 1):
        len_sum += result_clusters[i][1]
        if result_clusters[i + 1][1] / result_clusters[i][1] < len_ratio:
            main_groups = i + 1
            break
    return main_groups


def recommend_inflation(result_stat, nchrs, len_ratio):
    """
    根据各膨胀值对应的主要聚类数推荐最优膨胀参数。

    参数
    ----
    result_stat : list[tuple]
        [(inflation, main_groups), ...]
    nchrs : int
        预期染色体数
    len_ratio : float
        当前使用的长度比值阈值

    返回
    ----
    tuple[bool, float or None]
        (是否已给出建议, 推荐膨胀值)
    """
    separated = [
        (inf_, mg - nchrs)
        for inf_, mg in result_stat
        if mg >= int(nchrs)
    ]

    if separated:
        separated.sort(key=lambda x: x[0])
        rec = separated[0][0]
        logger.info('建议尝试膨胀值 {}（长度比值 = {}）'.format(rec, len_ratio))
        return True, rec
    else:
        if len_ratio > 0.5:
            logger.info('长度比值 {} 可能过严，尝试更低的值…'.format(len_ratio))
            return False, None
        else:
            logger.info(
                '部分染色体可能被合并分组（长度比值 = {}）。'
                '建议检查参数设置并尝试调整'.format(len_ratio)
            )
            return True, None


def run_mcl_clustering(
        link_matrix, bin_set, frag_len_dict, frag_index_dict,
        expansion, min_inflation, max_inflation, inflation_step,
        max_iter, pruning, fa_dict, nchrs, dense_matrix):
    """
    以不同膨胀值依次运行 Markov 聚类，并推荐最优膨胀参数。

    参数
    ----
    link_matrix : array or sparse matrix
        邻接矩阵
    bin_set : set
        所有 bin 名称集合
    frag_len_dict : dict
        {片段名: 长度}
    frag_index_dict : dict
        {片段名: 矩阵索引}
    expansion : int
        展开参数
    min_inflation, max_inflation : float
        膨胀值搜索范围
    inflation_step : float
        膨胀值步长
    max_iter : int
        最大迭代次数
    pruning : float
        剪枝阈值
    fa_dict : dict
        contig 信息字典
    nchrs : int
        预期染色体数
    dense_matrix : bool
        是否使用密集矩阵

    返回
    ----
    tuple
        (result_clusters_list, mcl_nrounds, rec_inf)
    """
    logger.info('执行 Markov 聚类…')

    rec_inf = 0
    idx_frag = {i: f for f, i in frag_index_dict.items()}

    start = Decimal(str(min_inflation))
    step  = Decimal(str(inflation_step))
    end   = Decimal(str(max_inflation)) + step

    # 预归一化并预展开
    matrix = normalize(link_matrix, norm='l1', axis=0)
    mem = get_matrix_memory_size(matrix)
    print('链接矩阵内存：{:.2f} MB'.format(mem / 1024 / 1024))

    if not dense_matrix:
        matrix = mkl_matrix_power(matrix, expansion)
    else:
        matrix = matrix_power(matrix, expansion)

    result_clusters_list = []
    mcl_nrounds = 0

    for inflation in arange(start, end, step):
        result_matrix = mcl(matrix, expansion, float(inflation), max_iter, pruning, dense_matrix)
        mcl_nrounds += 1

        clusters = interpret_result(result_matrix, dense_matrix)
        if not clusters:
            logger.info('膨胀值 {} 的结果存在片段缺失/重复，跳过输出'.format(inflation))
            continue

        result_clusters = defaultdict(lambda: [[], 0])
        ctg_clusters = defaultdict(dict)

        for n, indexes in enumerate(clusters):
            for i in indexes:
                frag = idx_frag[i]
                if frag in bin_set:
                    flen = frag_len_dict[frag]
                    ctg = frag.rsplit('_bin', 1)[0]
                    ctg_clusters[ctg][n] = ctg_clusters[ctg].get(n, 0) + flen
                else:
                    result_clusters[n][0].append(frag)
                    result_clusters[n][1] += fa_dict[frag][1]

        if ctg_clusters:
            for ctg, stats in ctg_clusters.items():
                best = max(stats, key=stats.get)
                result_clusters[best][0].append(ctg)
                result_clusters[best][1] += fa_dict[ctg][1]

        result_clusters = sorted(result_clusters.values(), key=lambda x: x[1], reverse=True)

        outdir = 'inflation_{}'.format(inflation)
        os.makedirs(outdir, exist_ok=True)

        with open('{0}/mcl_{0}.clusters.txt'.format(outdir), 'w') as fout:
            fout.write('#Group\tnContigs\tContigs\n')
            for n, (ctgs, glen) in enumerate(result_clusters, 1):
                result_clusters[n - 1][0].sort(key=lambda x: fa_dict[x][1], reverse=True)
                fout.write('group{}_{}bp\t{}\t{}\n'.format(
                    n, glen, len(ctgs), ' '.join(ctgs)))

        for n, (ctgs, glen) in enumerate(result_clusters, 1):
            with open('{}/group{}_{}bp.txt'.format(outdir, n, glen), 'w') as fout:
                fout.write('#Contig\tRECounts\tLength\n')
                for ctg in ctgs:
                    length, re_cnt = fa_dict[ctg][1:3]
                    fout.write('{}\t{}\t{}\n'.format(ctg, re_cnt, length))

        result_clusters_list.append((inflation, result_clusters))

    max_nclusters = max(len(rc) for _, rc in result_clusters_list)
    if max_nclusters < nchrs:
        logger.warning(
            '最大聚类数（{}）少于预期染色体数（{}），建议尝试更高的膨胀值'.format(
                max_nclusters, nchrs)
        )
    else:
        for lr in (0.75, 0.70, 0.65, 0.60, 0.55, 0.50, 0.40):
            stat = [(inf_, get_main_groups(rc, lr)) for inf_, rc in result_clusters_list]
            done, rec_inf = recommend_inflation(stat, nchrs, lr)
            if done:
                break

    return result_clusters_list, mcl_nrounds, rec_inf


# ─────────────────────────────────────────────────────────────
# 重分配辅助函数
# ─────────────────────────────────────────────────────────────

def add_ungrouped_ctgs(fa_dict, ctg_group_dict):
    """将未被分组的 contig 标记为 'ungrouped'。"""
    for ctg in fa_dict:
        if ctg not in ctg_group_dict:
            ctg_group_dict[ctg] = 'ungrouped'


def parse_link_dict(link_dict, ctg_group_dict):
    """
    统计每条 contig 与各聚类之间的链接总数。

    参数
    ----
    link_dict : dict
        完整的 contig 间链接字典
    ctg_group_dict : dict
        {contig: 所属聚类编号或 'ungrouped'}

    返回
    ----
    dict
        {contig: {聚类编号: 链接数}}
    """
    result = defaultdict(dict)

    def add(ctg, group, links):
        if group != 'ungrouped':
            result[ctg][group] = result[ctg].get(group, 0) + links

    for (ci, cj), links in link_dict.items():
        gi, gj = ctg_group_dict[ci], ctg_group_dict[cj]
        add(ci, gj, links)
        add(cj, gi, links)

    return result


def cal_link_density(max_group, current_group, max_links, group_RE_sites, ctg_RE_sites):
    """
    计算 contig 与目标聚类之间的链接密度。

    参数
    ----
    max_group : int
        链接数最多的聚类
    current_group : int or str
        contig 当前所在聚类
    max_links : int
        到目标聚类的链接数
    group_RE_sites : int
        目标聚类的 RE 位点总数
    ctg_RE_sites : int
        当前 contig 的 RE 位点数

    返回
    ----
    float
        链接密度
    """
    if max_group == current_group:
        return max_links / group_RE_sites
    return max_links / (group_RE_sites + ctg_RE_sites - 1)


# ─────────────────────────────────────────────────────────────
# 统计输出
# ─────────────────────────────────────────────────────────────

def output_statistics(fa_dict, link_dict, result_clusters_list):
    """
    生成各膨胀值结果的统计文件，用于后续重分配步骤的阈值选择。

    统计内容包括：RE 位点分布、Hi-C 链接数、链接密度、链接密度比值。
    若 matplotlib 可用，还会生成可视化图像。

    参数
    ----
    fa_dict : dict
        contig 信息字典
    link_dict : dict
        完整的 contig 间链接字典
    result_clusters_list : list
        [(inflation, result_clusters), ...]
    """

    def generate_axes(sorted_list):
        n_dict = OrderedDict({0: 0})
        l_dict = OrderedDict({0: 0})
        last = 0
        for ctg, val in sorted_list:
            if val in n_dict:
                n_dict[val] += 1
                l_dict[val] += fa_dict[ctg][1]
            else:
                n_dict[val] = n_dict[last] + 1
                l_dict[val] = l_dict[last] + fa_dict[ctg][1]
                last = val

        x, y1, y2 = [], [], []
        for k, v in n_dict.items():
            x.append(k)
            y1.append(v / total_n * 100)
            y2.append((total_len - l_dict[k]) / total_len * 100)
        return x, y1, y2

    def write_file(x, y1, y2, title, inflation):
        with open('inflation_{}/{}_statistics.txt'.format(inflation, title), 'w') as fh:
            fh.write('{}\tFiltered_ctg_n\tRest_ctg_len\n'.format(title))
            for i, v in enumerate(x):
                fh.write('>{}\t{}\t{}\n'.format(v, y1[i], y2[i]))

    logger.info('生成统计文件以辅助重分配步骤…')

    total_n = len(fa_dict)
    total_len = 0
    re_list = []
    for ctg, info in fa_dict.items():
        total_len += info[1]
        re_list.append((ctg, info[2]))
    re_list.sort(key=lambda x: x[1])
    xRE, y1RE, y2RE = generate_axes(re_list)

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import warnings
        warnings.filterwarnings('ignore', category=UserWarning)
        HAS_MPL = True
    except Exception:
        HAS_MPL = False
        logger.warning('matplotlib 未正确安装，将跳过可视化图像生成')

    for inflation, clusters in result_clusters_list:
        write_file(xRE, y1RE, y2RE, 'RE_site_threshold', inflation)

        ctg_grp = {}
        grp_re = {}
        for n, (ctgs, _) in enumerate(clusters):
            grp_re[n] = 1
            for c in ctgs:
                ctg_grp[c] = n
                grp_re[n] += fa_dict[c][2] - 1

        add_ungrouped_ctgs(fa_dict, ctg_grp)
        grp_links = parse_link_dict(link_dict, ctg_grp)

        link_list, density_list, ratio_list = [], [], []

        for ctg in fa_dict:
            if ctg in grp_links:
                sorted_gl = sorted(grp_links[ctg].items(), key=lambda x: x[1], reverse=True)
                mg, ml = sorted_gl[0]
                link_list.append((ctg, ml))
                cg = ctg_grp[ctg]
                md = cal_link_density(mg, cg, ml, grp_re[mg], fa_dict[ctg][2])
                density_list.append((ctg, md))
                if len(grp_re) > 1:
                    other_sum = sum(
                        cal_link_density(g, cg, lk, grp_re[g], fa_dict[ctg][2])
                        for g, lk in sorted_gl[1:]
                    )
                    avg_other = other_sum / (len(grp_re) - 1)
                else:
                    avg_other = 0
                ratio_list.append((ctg, md / avg_other if avg_other else 1_000_000))
            else:
                link_list.append((ctg, 0))
                density_list.append((ctg, 0))
                ratio_list.append((ctg, 0))

        link_list.sort(key=lambda x: x[1])
        xL, y1L, y2L = generate_axes(link_list)
        write_file(xL, y1L, y2L, 'Link_threshold', inflation)

        density_list.sort(key=lambda x: x[1])
        xD, y1D, y2D = generate_axes(density_list)
        write_file(xD, y1D, y2D, 'Link_density_threshold', inflation)

        ratio_list.sort(key=lambda x: x[1])
        xR, y1R, y2R = generate_axes(ratio_list)
        write_file(xR, y1R, y2R, 'Link_density_ratio_threshold', inflation)

        if HAS_MPL:
            fig, axes = plt.subplots(2, 2, figsize=(8, 7))

            def _plot(ax, x, y1, y2, xlim, ylim_b, title, xlabel):
                ax.plot(x, y1, 'b')
                ax.tick_params(axis='y', colors='b')
                ax.set_xlim([0, xlim])
                ax.set_ylim([0, ylim_b])
                ax.set_ylabel('过滤 contig 比例 (%)', color='b')
                ax.set_title(title)
                ax.set_xlabel(xlabel)
                ax2 = ax.twinx()
                ax2.plot(x, y2, 'r')
                ax2.tick_params(axis='y', colors='r')
                ax2.set_ylim([90, 100])
                ax2.set_ylabel('剩余 contig 长度 (%)', color='r')

            _plot(axes[0][0], xRE, y1RE, y2RE, 500, 50, 'RE 位点阈值', 'RE 位点数')
            _plot(axes[0][1], xL,  y1L,  y2L,  500, 50, 'Hi-C 链接阈值', '最优组链接数')
            _plot(axes[1][0], xD,  y1D,  y2D,  0.001, 50, '链接密度阈值', '最优组链接密度')
            _plot(axes[1][1], xR,  y1R,  y2R,  20,  50, '链接密度比值阈值', '密度比值（最优/均值）')

            fig.tight_layout(w_pad=1, h_pad=1)
            plt.savefig('inflation_{}/statistics.pdf'.format(inflation))
            plt.close()


# ─────────────────────────────────────────────────────────────
# 参数校验
# ─────────────────────────────────────────────────────────────

def check_param(param, string, suffix, true_suffix=''):
    """
    校验字符串型参数的合法性，支持纯数值和带后缀（如 'X'）两种模式。

    参数
    ----
    param : str
        参数名称（用于错误提示）
    string : str
        参数值字符串
    suffix : set
        合法的后缀字符集合
    true_suffix : str
        递归时传递已识别的后缀

    返回
    ----
    tuple
        (数值, 后缀字符串)

    异常
    ----
    RuntimeError
        参数不合法时抛出
    """
    if len(string) > 1:
        if suffix and string[-1] in suffix:
            return check_param(param, string[:-1], None, string[-1])
        try:
            num = float(string)
            if not true_suffix and not 0 <= num <= 1:
                logger.error('参数 {} {} 非法'.format(param, string + true_suffix))
            else:
                return num, true_suffix
        except ValueError:
            logger.error('参数 {} {} 非法'.format(param, string + true_suffix))
    elif len(string) == 1:
        try:
            num = float(string)
            if not true_suffix and not 0 <= num <= 1:
                logger.error('参数 {} {} 非法'.format(param, string + true_suffix))
            else:
                return num, true_suffix
        except ValueError:
            logger.error('参数 {} {} 非法'.format(param, string + true_suffix))
    else:
        logger.error('参数 {} 为空'.format(param))

    raise RuntimeError('参数校验失败')


# ─────────────────────────────────────────────────────────────
# 格式检测与内存统计
# ─────────────────────────────────────────────────────────────

def detect_format(args):
    """根据文件扩展名自动判断 Hi-C 比对文件的格式。"""
    if args.alignments.endswith('.bam'):
        args.aln_format = 'bam'
        logger.info('检测到 BAM 格式的比对文件')
    elif args.alignments.endswith('.pairs'):
        args.aln_format = 'pairs'
        logger.info('检测到 pairs 格式的比对文件')
    elif args.alignments.endswith('.pairs.gz'):
        args.aln_format = 'bgzipped_pairs'
        logger.info('检测到 bgzipped pairs 格式的比对文件')
    else:
        raise RuntimeError('无法识别的 Hi-C 比对文件格式')


def get_matrix_memory_size(matrix):
    """
    估算矩阵对象的内存占用量（字节）。

    对稀疏矩阵，累加 data、indices、indptr 各自的内存；
    对密集矩阵，直接调用 sys.getsizeof。

    参数
    ----
    matrix : array or sparse matrix
        待估算的矩阵

    返回
    ----
    int
        内存占用（字节）
    """
    if sp.issparse(matrix):
        return (
            sys.getsizeof(matrix)
            + sys.getsizeof(matrix.data)
            + sys.getsizeof(matrix.indices)
            + sys.getsizeof(matrix.indptr)
        )
    return sys.getsizeof(matrix)


# ─────────────────────────────────────────────────────────────
# 参数解析
# ─────────────────────────────────────────────────────────────

def parse_arguments():
    """解析命令行参数并返回 Namespace 对象。"""

    parser = argparse.ArgumentParser(prog='conhic cluster')

    # 输入文件与流程控制
    inp = parser.add_argument_group('>>> 输入文件与流程控制参数')
    inp.add_argument('fasta', help='草图基因组 FASTA 文件')
    inp.add_argument('alignments', help='经过过滤的 Hi-C 比对文件（BAM 或 pairs 格式，请勿按坐标排序）')
    inp.add_argument('nchrs', type=int, help='预期染色体数目')
    inp.add_argument(
        '--aln_format', choices={'bam', 'pairs', 'bgzipped_pairs', 'auto'},
        default='auto', help='Hi-C 比对文件格式，默认：%(default)s')
    inp.add_argument(
        '--RE', default='GATC',
        help='限制酶识别位点，多位点用逗号分隔，默认：%(default)s')

    # 片段过滤参数
    flt = parser.add_argument_group('>>> 聚类前预处理（片段与链接过滤）参数')
    flt.add_argument(
        '--Nx', type=int, default=80,
        help='保留长度满足 Nx 的 contig，默认：%(default)s')
    flt.add_argument(
        '--RE_site_cutoff', type=int, default=5,
        help='最低 RE 位点数阈值，默认：%(default)s')
    flt.add_argument(
        '--density_lower', default='0.2X',
        help='链接密度下限，支持分数模式（0~1）和倍数模式（以 X 结尾），默认：%(default)s')
    flt.add_argument(
        '--density_upper', default='1.9X',
        help='链接密度上限，默认：%(default)s')
    flt.add_argument(
        '--read_depth_upper', default='1.5X',
        help='测序深度上限，默认：%(default)s')
    flt.add_argument(
        '--topN', type=int, default=10,
        help='秩和计算时使用的最近邻片段数，默认：%(default)s')
    flt.add_argument(
        '--rank_sum_hard_cutoff', type=int, default=0,
        help='秩和硬过滤阈值，默认：%(default)s（禁用）')
    flt.add_argument(
        '--rank_sum_upper', default='1.5X',
        help='秩和上限，默认：%(default)s')
    flt.add_argument(
        '--remove_allelic_links', type=int, default=0,
        help='识别并移除等位 contig 之间的 Hi-C 链接；值应为倍性数（≥2），默认禁用')
    flt.add_argument(
        '--concordance_ratio_cutoff', type=float, default=0.2,
        help='等位 contig 识别的一致性比率阈值，默认：%(default)s')
    flt.add_argument(
        '--nwindows', type=int, default=50,
        help='一致性比率计算中的分窗数目，默认：%(default)s')
    flt.add_argument(
        '--remove_concentrated_links', default=False, action='store_true',
        help='移除高度集中于局部区域的 Hi-C 链接，默认：%(default)s')
    flt.add_argument(
        '--max_read_pairs', type=int, default=200,
        help='等位 contig 识别的最大 read 对数，默认：%(default)s')
    flt.add_argument(
        '--min_read_pairs', type=int, default=20,
        help='等位 contig 识别的最小 read 对数，默认：%(default)s')
    flt.add_argument(
        '--phasing_weight', type=float, default=1.0,
        help='单倍型信息权重，默认：%(default)s（1.0 表示完全移除跨单倍型链接）')

    # 矩阵构建与 Markov 聚类参数
    mcl = parser.add_argument_group('>>> 邻接矩阵构建与 Markov 聚类参数')
    mcl.add_argument(
        '--bin_size', type=int, default=-1,
        help='分箱大小（kbp），超过此长度的 contig 将被切割；-1 表示自动计算，0 表示禁用，默认：%(default)s')
    mcl.add_argument(
        '--flank', type=int, default=500,
        help='仅使用 contig 两端侧翼区域的链接构建矩阵，单位 kbp，默认：%(default)s')
    mcl.add_argument(
        '--normalize_by_nlinks', default=False, action='store_true',
        help='按总链接数归一化侧翼链接，默认：%(default)s')
    mcl.add_argument(
        '--expansion', type=int, default=2,
        help='Markov 聚类展开参数，默认：%(default)s')
    mcl.add_argument(
        '--min_inflation', type=float, default=1.2,
        help='膨胀参数搜索下限，默认：%(default)s')
    mcl.add_argument(
        '--max_inflation', type=float, default=1.8,
        help='膨胀参数搜索上限，默认：%(default)s')
    mcl.add_argument(
        '--inflation_step', type=float, default=0.1,
        help='膨胀参数步长，默认：%(default)s')
    mcl.add_argument(
        '--max_iter', type=int, default=200,
        help='每次 Markov 聚类的最大迭代次数，默认：%(default)s')
    mcl.add_argument(
        '--pruning', type=float, default=0.0001,
        help='剪枝阈值，默认：%(default)s')
    mcl.add_argument(
        '--skip_clustering', default=False, action='store_true',
        help='跳过 Markov 聚类步骤，默认：%(default)s')

    # 性能参数
    perf = parser.add_argument_group('>>> 性能参数')
    perf.add_argument(
        '--threads', type=int, default=8,
        help='读取 BAM 文件的线程数，默认：%(default)s')
    perf.add_argument(
        '--dense_matrix', default=False, action='store_true',
        help='使用密集矩阵模式（默认使用稀疏矩阵 + MKL 加速），默认：%(default)s')

    # 日志参数
    log = parser.add_argument_group('>>> 日志参数')
    log.add_argument(
        '--verbose', default=False, action='store_true',
        help='输出详细日志，默认：%(default)s')

    args = parser.parse_args()

    check_param('--density_lower',   args.density_lower,   {'X', 'x'})
    check_param('--density_upper',   args.density_upper,   {'X', 'x'})
    check_param('--read_depth_upper', args.read_depth_upper, {'X', 'x'})
    check_param('--rank_sum_upper',  args.rank_sum_upper,  {'X', 'x'})

    return args


# ─────────────────────────────────────────────────────────────
# 主流程
# ─────────────────────────────────────────────────────────────

def run(args, log_file=None, random_times=10, origin_path=None):
    """
    ConHiC 聚类主流程。

    依次执行以下步骤：
    1. 读取并统计 FASTA 序列信息；
    2. 对多种 bin_size 分别构建邻接矩阵；
    3. 解析 Hi-C 比对，过滤片段；
    4. 运行 Markov 聚类并推荐最优膨胀参数。

    参数
    ----
    args : argparse.Namespace
        命令行参数
    log_file : str, optional
        附加日志文件路径
    random_times : int
        随机采样组数（用于小 bin_size 的稳定性评估）
    origin_path : str, optional
        若已有中间结果，可直接从此路径加载（跳过重新解析）
    """
    if log_file:
        fh = logging.FileHandler(log_file, 'w')
        fh.setFormatter(logging.Formatter(
            fmt='%(asctime)s <%(filename)s> [%(funcName)s] %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        ))
        logger.addHandler(fh)

    start_time = time.time()
    logger.info('ConHiC 启动，版本：{} （更新日期：{}}）'.format(__version__, __update_time__))
    logger.info('Python 版本：{}'.format(sys.version.replace('\n', '')))
    logger.info('命令：{}'.format(' '.join(sys.argv)))

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    if not INTEL_MKL:
        logger.warning('sparse_dot_mkl 或 Intel MKL 未正确安装，将以密集矩阵模式运行')
        args.dense_matrix = True
    elif args.dense_matrix:
        logger.warning('已指定 --dense_matrix，将以密集矩阵模式运行')

    if args.aln_format == 'auto':
        detect_format(args)

    # 读取基因组 FASTA
    fa_dict = parse_fasta(args.fasta, RE=args.RE)
    pos_int_type, dist_int_type = determine_int_type(fa_dict)

    import numpy as np
    import matplotlib.pyplot as plt

    print(len(fa_dict), ' 条 contig')

    lengths = [info[1] for info in fa_dict.values()]
    mean_len = np.mean(lengths)

    # 构造多组 bin_size 进行多尺度评估
    bin_size_list = [int(mean_len * m) for m in (20, 19, 18, 17, 16)]

    read_depth_dict = {}
    whitelist = set()
    args.whitelist = whitelist
    rec_infs = {}

    for b_s in bin_size_list:
        os.makedirs(str(b_s), exist_ok=True)
        os.chdir(str(b_s))

        fa_dict = parse_fasta(args.fasta, RE=args.RE)
        (_, bin_set, bin_size, frag_len_dict, Nx_frag_set,
         RE_site_dict, split_ctg_set, fa_bin_dict, fa_bin_num) = stat_fragments(
            fa_dict, args.RE, read_depth_dict, whitelist,
            nchrs=args.nchrs, flank=args.flank, Nx=args.Nx, bin_size=b_s
        )

        output_pickle(fa_bin_dict, 'fa_bin_dict', 'fa_bin_dict.pkl')
        output_pickle(fa_bin_num, 'fa_bin_num', 'fa_bin_num.pkl')

        all_results = generate_multiple_results(
            fa_bin_dict, select_ratio=0.25, n_groups=random_times, base_seed=42
        )
        print('随机采样 bin 数：', len(all_results[0]))

        vals = list(fa_bin_num.values())
        plt.hist(vals, bins=10, color='blue', edgecolor='black')
        plt.xlabel('Bin 数量')
        plt.ylabel('频次')
        plt.title('各 contig 的 bin 数量分布')
        plt.savefig('histogram2.png', dpi=300)
        plt.close()

        # 根据比对文件格式选取对应生成器
        if split_ctg_set:
            if args.aln_format == 'bam':
                format_options = [b'filter=flag.read1']
                alignments = bam_generator(args.alignments, args.threads, format_options)
            else:
                alignments = pairs_generator(args.alignments, args.aln_format)
        else:
            if args.aln_format == 'bam':
                format_options = [b'filter=flag.read1 && refid != mrefid']
                alignments = bam_generator(args.alignments, args.threads, format_options)
            else:
                alignments = pairs_generator_inter_ctgs(args.alignments, args.aln_format)

        # 解析比对文件并构建链接字典
        if True:
            if split_ctg_set:
                (full_link_dict, flank_link_dict, HT_link_dict,
                 clm_dict, frag_link_dict, ctg_coord_dict, ctg_pair_to_frag) = parse_alignments(
                    alignments, fa_dict, args, bin_size, frag_len_dict,
                    Nx_frag_set, split_ctg_set, pos_int_type, dist_int_type
                )
            else:
                (full_link_dict, flank_link_dict, HT_link_dict,
                 clm_dict, frag_link_dict, ctg_coord_dict) = parse_alignments_for_ctgs(
                    alignments, fa_dict, args, frag_len_dict, Nx_frag_set, pos_int_type, dist_int_type
                )

            output_pickle(flank_link_dict, 'flank_link_dict', 'flank_link_dict.pkl')
            output_pickle(frag_link_dict,  'frag_link_dict',  'frag_link_dict.pkl')
            output_pickle(HT_link_dict,    'HT_link_dict',    'HT_links.pkl')
            del HT_link_dict
            gc.collect()

            output_clm(clm_dict)
            del clm_dict
            gc.collect()

            # 移除集中链接对 full_link_dict 的影响
            if args.remove_concentrated_links:
                for pair, data in ctg_coord_dict.items():
                    if isinstance(data, list):
                        full_link_dict[pair] *= data[1]

            # 片段过滤
            filtered_frags = filter_fragments(
                Nx_frag_set, RE_site_dict, args.RE_site_cutoff, frag_link_dict,
                args.density_lower, args.density_upper, args.topN, args.rank_sum_upper,
                args.rank_sum_hard_cutoff, flank_link_dict, read_depth_dict,
                args.read_depth_upper, whitelist
            )
            output_pickle(filtered_frags, 'filtered_frags', 'filtered_frags.pkl')

            # 移除等位链接
            if args.remove_allelic_links:
                if split_ctg_set:
                    filtered_frags = remove_allelic_HiC_links(
                        fa_dict, ctg_coord_dict, full_link_dict, args,
                        flank_link_dict, filtered_frags, ctg_pair_to_frag
                    )
                else:
                    filtered_frags = remove_allelic_HiC_links(
                        fa_dict, ctg_coord_dict, full_link_dict, args,
                        flank_link_dict, filtered_frags
                    )
                del ctg_coord_dict
                gc.collect()

            output_pickle(full_link_dict, 'full_link_dict', 'full_links.pkl')

        mem = get_matrix_memory_size(full_link_dict)
        print('full_link_dict 内存：{:.2f} MB'.format(mem / 1024 / 1024))

        if args.normalize_by_nlinks:
            normalize_by_nlinks(flank_link_dict, frag_link_dict)
        del frag_link_dict

        # 对于较小的 bin_size，使用随机采样策略评估稳定性
        if b_s <= int(mean_len / 6):
            rec_infs[str(b_s)] = []

            for i, sel_bins in enumerate(all_results):
                os.makedirs(str(i), exist_ok=True)
                os.chdir(str(i))

                filtered_dict = remove_selected_bins_from_flank_links(flank_link_dict, sel_bins)
                mat, frag_idx = dict_to_matrix(
                    filtered_dict, filtered_frags,
                    dense_matrix=args.dense_matrix, add_self_loops=True
                )
                del filtered_dict
                gc.collect()

                logger.info('链接矩阵构建耗时：{:.1f}s'.format(time.time() - start_time))

                if not args.skip_clustering:
                    clusters_list, nrounds, rec_inf = run_mcl_clustering(
                        mat, bin_set, frag_len_dict, frag_idx,
                        args.expansion, args.min_inflation, args.max_inflation,
                        args.inflation_step, args.max_iter, args.pruning,
                        fa_dict, args.nchrs, args.dense_matrix
                    )
                    rec_infs[str(b_s)].append(rec_inf)
                    output_statistics(fa_dict, full_link_dict, clusters_list)

                logger.info('完成，总耗时：{:.1f}s'.format(time.time() - start_time))
                os.chdir('..')

            # 取众数作为该 bin_size 的最终推荐值
            vals, cnts = np.unique(rec_infs[str(b_s)], return_counts=True)
            rec_infs[str(b_s)] = vals[np.argmax(cnts)]

        else:
            # 对于较大的 bin_size，直接运行聚类
            mat, frag_idx = dict_to_matrix(
                flank_link_dict, filtered_frags,
                dense_matrix=args.dense_matrix, add_self_loops=True
            )
            del filtered_frags, flank_link_dict
            gc.collect()

            logger.info('链接矩阵构建耗时：{:.1f}s'.format(time.time() - start_time))

            if not args.skip_clustering:
                clusters_list, nrounds, rec_inf = run_mcl_clustering(
                    mat, bin_set, frag_len_dict, frag_idx,
                    args.expansion, args.min_inflation, args.max_inflation,
                    args.inflation_step, args.max_iter, args.pruning,
                    fa_dict, args.nchrs, args.dense_matrix
                )
                rec_infs[str(b_s)] = rec_inf
                output_statistics(fa_dict, full_link_dict, clusters_list)

            logger.info('完成，总耗时：{:.1f}s'.format(time.time() - start_time))

        os.chdir('..')

    logger.info('各 bin_size 推荐膨胀值：{}'.format(rec_infs))
    return rec_infs, bin_size_list


def main():
    """程序入口：解析参数并启动主流程。"""
    args = parse_arguments()
    run(args, 'ConHiC_cluster.log')


if __name__ == '__main__':
    main()

