import matplotlib.pyplot as plt
import numpy as np
import os
import sys
from matplotlib.ticker import ScalarFormatter, FuncFormatter

# 定义文件路径和对应的标签
file_paths = {
    'AP': '/home/chenzh/HapHiC/data/ap/shuffled_ap_split_genome.fa',
    'C88_2M': '/home/chenzh/HapHiC/data/C88/sim/C88.v1.split_2000000_0.2.fa',
    'C88_0.8M': '/home/chenzh/HapHiC/data/C88/sim2/C88.v1.split_800000_0.2.fa',
    'CS': '/home/chenzh/HapHiC/data/CS/shuffled_cs_split_genome.fasta',
    'SS': '/home/chenzh/HapHiC/data/Saccharum_spontaneum_Np-X/shuffled_Np-X_split_genome.fa',
    'XJDY': '/home/chenzh/HapHiC/data/XinJiangDaYe/shuffled_xjdy_split_genome.fa',
    'ZM-4': '/home/chenzh/HapHiC/data/ZM-4/shuffled_zm-4_split_genome.fa'
}

def read_fasta_lengths(file_path):
    """读取FASTA文件并返回所有序列的长度列表"""
    lengths = []
    current_length = 0
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):  # 序列头
                if current_length > 0:  # 保存上一个序列的长度
                    lengths.append(current_length)
                current_length = 0  # 重置计数器
            else:  # 序列行
                current_length += len(line)
        # 保存最后一个序列的长度
        if current_length > 0:
            lengths.append(current_length)
    return lengths

def sci_notation_millions(x, pos):
    """将坐标值格式化为10^6的科学计数法"""
    return f'{x/1e6:.1f}'

# 创建图形和子图 - 改为2行4列
fig, axes = plt.subplots(2, 4, figsize=(20, 10))
axes = axes.flatten()

# 设置全局参数
plt.rcParams['font.size'] = 10
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['axes.labelsize'] = 10

# 用于存储所有数据的统计信息
all_stats = []

# 对每个文件进行处理
for idx, (sample_name, file_path) in enumerate(file_paths.items()):
    if idx >= len(axes):
        break
    
    ax = axes[idx]
    
    # 检查文件是否存在
    if not os.path.exists(file_path):
        print(f"警告: 文件不存在 - {file_path}")
        ax.text(0.5, 0.5, f"File not found\n{os.path.basename(file_path)}", 
                ha='center', va='center', transform=ax.transAxes)
        continue
    
    # 读取contig长度
    try:
        lengths = read_fasta_lengths(file_path)
        if not lengths:
            ax.text(0.5, 0.5, "No sequences found", ha='center', va='center', transform=ax.transAxes)
            continue
            
        lengths = np.array(lengths)
        print(f"{sample_name}: 读取到 {len(lengths)} 个contigs")
        
    except Exception as e:
        print(f"读取 {sample_name} 时出错: {e}")
        ax.text(0.5, 0.5, f"Error reading file\n{os.path.basename(file_path)}", 
                ha='center', va='center', transform=ax.transAxes)
        continue
    
    # 计算统计信息
    mean_len = np.mean(lengths)
    median_len = np.median(lengths)
    min_len = np.min(lengths)
    max_len = np.max(lengths)
    total_contigs = len(lengths)
    below_mean_ratio = np.sum(lengths < mean_len) / total_contigs
    
    # 存储统计信息
    all_stats.append({
        'sample': sample_name,
        'mean': mean_len,
        'median': median_len,
        'min': min_len,
        'max': max_len,
        'total': total_contigs,
        'below_mean_ratio': below_mean_ratio
    })
    
    # 绘制直方图
    n_bins = min(50, int(np.sqrt(len(lengths))))
    n_bins = max(n_bins, 10)
    
    ax.hist(lengths, bins=n_bins, alpha=0.7, color='steelblue', edgecolor='black')
    
    # 添加均值线
    ax.axvline(x=mean_len, color='red', linestyle='--', linewidth=1.5, 
               label=f'Mean: {mean_len:.0f} bp')
    
    # 添加标题和标签
    ax.set_title(f'{sample_name} (N={total_contigs:,})', fontweight='bold')
    ax.set_xlabel('Contig Length (bp)')
    ax.set_ylabel('Frequency')
    
    # 统一设置x轴为10^6科学计数法
    ax.xaxis.set_major_formatter(FuncFormatter(sci_notation_millions))
    ax.set_xlabel('Contig Length (×10⁶ bp)')  # 更新x轴标签
    
    # 设置x轴范围的下限为0
    ax.set_xlim(left=0)
    
    # 添加统计信息文本
    stats_text = f'Mean: {mean_len:,.0f} bp\nMedian: {median_len:,.0f} bp\n'
    stats_text += f'Min: {min_len:,.0f} bp\nMax: {max_len:,.0f} bp\n'
    stats_text += f'<Mean: {below_mean_ratio:.1%}'
    
    # 将文本放在合适的位置
    ax.text(0.98, 0.98, stats_text, transform=ax.transAxes,
            verticalalignment='top', horizontalalignment='right',
            fontsize=9, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    # 添加网格
    ax.grid(True, alpha=0.3, linestyle='--')
    
    # 添加图例
    ax.legend(loc='upper left')

# 隐藏最后一个空的子图（如果有的话）
if len(file_paths) < len(axes):
    for i in range(len(file_paths), len(axes)):
        axes[i].set_visible(False)

# 调整布局
plt.tight_layout()

# 添加总标题
fig.suptitle('Contig Length Distribution Across Samples', fontsize=16, fontweight='bold', y=1.02)

# 保存图形
output_file = 'contig_length_distributions_2x4.png'
plt.savefig(output_file, dpi=300, bbox_inches='tight')
print(f"\n图形已保存为: {output_file}")

# 显示图形
plt.show()

# 打印汇总统计
print("\n" + "="*80)
print("SUMMARY STATISTICS:")
print("="*80)
print(f"{'Sample':<10} {'Contigs':<10} {'Mean(bp)':<12} {'Median(bp)':<12} {'Min(bp)':<12} {'Max(bp)':<12} {'<Mean%':<10}")
print("-"*80)

for stats in all_stats:
    print(f"{stats['sample']:<10} {stats['total']:<10,} {stats['mean']:<12,.0f} "
          f"{stats['median']:<12,.0f} {stats['min']:<12,.0f} {stats['max']:<12,.0f} "
          f"{stats['below_mean_ratio']:<10.1%}")

print("="*80)