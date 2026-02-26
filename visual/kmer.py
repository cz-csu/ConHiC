import os
import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from collections import defaultdict
import matplotlib.gridspec as gridspec

def parse_chromosome(chrom_name):
    """从染色体名称中提取主染色体号（如从chr1_1提取chr1）"""
    match = re.match(r'(chr\d+)', chrom_name)
    if match:
        return match.group(1)
    return chrom_name

def get_dominant_chromosome(chrom_counts):
    """确定该group主要属于哪个染色体"""
    # 统计每个主染色体的总计数
    main_chrom_counts = defaultdict(int)
    total_count = 0
    
    for chrom, count in chrom_counts:
        main_chrom = parse_chromosome(chrom)
        if main_chrom.startswith('chr'):
            main_chrom_counts[main_chrom] += count
        total_count += count
    
    if not main_chrom_counts or total_count == 0:
        return "unreliable"
    
    # 找到计数最多的主染色体
    dominant_chrom = max(main_chrom_counts.items(), key=lambda x: x[1])[0]
    
    # 计算占比
    dominant_ratio = main_chrom_counts[dominant_chrom] / total_count
    
    # 如果占比太低，标记为unreliable
    if dominant_ratio < 0.5:
        return "unreliable"
    
    return dominant_chrom

def parse_file(filepath):
    """解析单个group文件"""
    data = []
    chrom_counts_per_interval = []
    
    with open(filepath, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 6:
                continue
                
            # 提取前5列
            group = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            hap = parts[3]
            score = float(parts[4])
            
            # 解析染色体计数（第6列）
            chrom_counts_str = parts[5].replace('dict_items([(', '').replace(')])', '')
            chrom_counts = []
            
            # 解析每个染色体计数对
            items = chrom_counts_str.split('), (')
            for item in items:
                item = item.strip("'()")
                if item:
                    sub_items = item.split("', ")
                    if len(sub_items) == 2:
                        chrom_name = sub_items[0].strip("'")
                        count = int(sub_items[1].strip("'"))
                        chrom_counts.append((chrom_name, count))
            
            data.append({
                'group': group,
                'start': start,
                'end': end,
                'hap': hap,
                'score': score,
                'chrom_counts': chrom_counts
            })
            chrom_counts_per_interval.append(chrom_counts)
    
    return data, chrom_counts_per_interval

def process_all_groups(data_dir):
    """处理所有group文件"""
    groups_data = {}
    group_dominant_chrom = {}
    
    # 获取所有group文件
    group_files = [f for f in os.listdir(data_dir) if f.startswith('group') and f.endswith('.txt')]
    
    for filename in group_files:
        group_name = filename.split('_')[0]
        filepath = os.path.join(data_dir, filename)
        
        # 解析文件
        data, chrom_counts_per_interval = parse_file(filepath)
        
        if not data:
            continue
            
        # 确定该group的主要染色体
        all_chrom_counts = []
        for chrom_counts in chrom_counts_per_interval:
            all_chrom_counts.extend(chrom_counts)
        
        dominant_chrom = get_dominant_chromosome(all_chrom_counts)
        
        groups_data[group_name] = data
        group_dominant_chrom[group_name] = dominant_chrom
    
    return groups_data, group_dominant_chrom

def create_chromosome_plot_data(groups_data, group_dominant_chrom):
    """为所有染色体创建绘图数据"""
    # 染色体顺序
    chromosomes = [f'chr{i}' for i in range(1, 13)]
    
    # 存储每个染色体的数据
    chrom_data = {}
    
    for chrom in chromosomes:
        # 获取所有属于该染色体的groups
        target_groups = [g for g, c in group_dominant_chrom.items() 
                        if c == chrom and g in groups_data]
        
        if not target_groups:
            continue
        
        # 排序groups
        target_groups.sort(key=lambda x: int(x.replace('group', '')))
        
        # 存储每个group的数据（使用所有区间）
        group_data_dict = {}
        
        for group in target_groups:
            data = groups_data[group]
            actual_intervals = len(data)
            
            # 创建数据数组 - 只需要hap信息，不需要score
            haps = np.zeros(actual_intervals, dtype=int)
            
            # 填充数据
            for j, interval in enumerate(data):
                # 将hap转换为数字编码
                hap = interval['hap']
                if hap == 'hap1':
                    haps[j] = 0  # 红色
                elif hap == 'hap2':
                    haps[j] = 1  # 蓝色
                elif hap == 'hap3':
                    haps[j] = 2  # 绿色
                elif hap == 'hap4':
                    haps[j] = 3  # 紫色
                elif 'chr' in hap and hap.startswith('chr'):
                    haps[j] = 4  # 棕色 (其他染色体)
                else:
                    haps[j] = 5  # 黄色 (不可靠)
            
            group_data_dict[group] = {
                'haps': haps,
                'num_intervals': actual_intervals
            }
        
        chrom_data[chrom] = {
            'groups': target_groups,
            'group_data': group_data_dict
        }
    
    return chrom_data

def get_color_for_hap(hap_type):
    """根据hap类型获取颜色（固定颜色，不使用score调整深浅）"""
    # 固定颜色
    base_colors = {
        0: (1.0, 0.0, 0.0),    # 红色 - hap1
        1: (0.0, 0.0, 1.0),    # 蓝色 - hap2
        2: (0.0, 1.0, 0.0),    # 绿色 - hap3
        3: (0.5, 0.0, 0.5),    # 紫色 - hap4
        4: (0.6, 0.4, 0.2),    # 棕色 - other_chrom
        5: (1.0, 1.0, 0.0)     # 黄色 - unreliable
    }
    
    return base_colors[hap_type]

def print_chromosome_groups_bash_format(group_dominant_chrom):
    """以bash数组格式输出每个染色体对应的group"""
    print("\n" + "="*60)
    print("染色体分组信息 (Bash数组格式):")
    print("="*60)
    
    # 按染色体号排序的groups字典
    chrom_groups_dict = {}
    
    # 初始化所有染色体（1-12）
    for i in range(1, 13):
        chrom_groups_dict[str(i)] = []
    
    # 收集每个染色体的group编号
    for group, chrom in group_dominant_chrom.items():
        if chrom.startswith('chr'):
            # 提取染色体编号
            chrom_num = chrom.replace('chr', '')
            if chrom_num.isdigit() and 1 <= int(chrom_num) <= 12:
                # 提取group编号（去掉"group"前缀）
                group_num = group.replace('group', '')
                if group_num.isdigit():
                    chrom_groups_dict[chrom_num].append(group_num)
    
    # 按group编号排序每个染色体内的groups
    for chrom in chrom_groups_dict:
        chrom_groups_dict[chrom].sort(key=int)
    
    # 打印bash数组格式
    print("declare -A chr_groups=(")
    for chrom_num in sorted(chrom_groups_dict.keys(), key=int):
        groups_list = chrom_groups_dict[chrom_num]
        if groups_list:  # 只输出有groups的染色体
            groups_str = "|".join(groups_list)
            print(f"    [{chrom_num}]=\"{groups_str}\"")
    print(")")
    
    # 另外打印一个更易读的格式
    print("\n" + "-"*60)
    print("易读格式:")
    print("-"*60)
    for chrom_num in sorted(chrom_groups_dict.keys(), key=int):
        groups_list = chrom_groups_dict[chrom_num]
        if groups_list:
            groups_str = "|".join(groups_list)
            print(f"chr{chrom_num}: {groups_str}")
    
    # 输出未分配的groups（如果有）
    unassigned_groups = []
    for group, chrom in group_dominant_chrom.items():
        if chrom == "unreliable" or not chrom.startswith('chr'):
            unassigned_groups.append(group.replace('group', ''))
    
    if unassigned_groups:
        print("\n" + "-"*60)
        print(f"未分配的groups ({len(unassigned_groups)}个):")
        print("-"*60)
        unassigned_groups.sort(key=int)
        print("|".join(unassigned_groups))
    
    # 统计信息
    total_groups = 0
    for groups_list in chrom_groups_dict.values():
        total_groups += len(groups_list)
    
    print("\n" + "="*60)
    print(f"统计信息:")
    print(f"- 总group数: {total_groups}")
    for chrom_num in sorted(chrom_groups_dict.keys(), key=int):
        groups_list = chrom_groups_dict[chrom_num]
        if groups_list:
            print(f"- chr{chrom_num}: {len(groups_list)} groups")
    print("="*60)

def plot_all_chromosomes_combined(groups_data, group_dominant_chrom, output_dir='combined_heatmaps'):
    """将所有染色体的热力图合并到一张图上"""
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    # 创建绘图数据
    chrom_data = create_chromosome_plot_data(groups_data, group_dominant_chrom)
    
    if not chrom_data:
        print("没有找到有效的数据！")
        return
    
    # 染色体顺序
    chromosomes = [f'chr{i}' for i in range(1, 13)]
    
    # 创建图形 - 增加高度以容纳下移的chr7-chr12
    fig = plt.figure(figsize=(30, 24))
    
    # 创建网格布局
    n_rows = 6  # 逻辑行数
    n_cols = 6  # 每行6个染色体
    
    # ================= 修改处 1 =================
    # 调整 height_ratios：
    # [1, 1] -> 第0-1行 (Chr1-6 使用)
    # [0.15] -> 第2行 (作为间隔 Gap)。将其从 1 改为 0.15，大幅减小中间空隙，使下方图表上移。
    # [1, 1] -> 第3-4行 (Chr7-12 使用)
    # [1]    -> 第5行 (底部余量)
    gs = gridspec.GridSpec(n_rows, n_cols, figure=fig, height_ratios=[1, 1, 0.15, 1, 1, 1])
    # ============================================
    
    # 颜色定义和标签
    hap_colors = ['red', 'blue', 'green', 'purple', 'brown', 'yellow']
    hap_rgb = [
        (1.0, 0.0, 0.0),    # 红色
        (0.0, 0.0, 1.0),    # 蓝色
        (0.0, 1.0, 0.0),    # 绿色
        (0.5, 0.0, 0.5),    # 紫色
        (0.6, 0.4, 0.2),    # 棕色
        (1.0, 1.0, 0.0)     # 黄色
    ]
    hap_labels = ['hap1', 'hap2', 'hap3', 'hap4', 'other_chrom', 'unreliable']
    
    # 计算所有染色体中最大的group数量（用于统一y轴范围）
    max_groups_per_chrom = 0
    for chrom in chromosomes:
        if chrom in chrom_data:
            num_groups = len(chrom_data[chrom]['groups'])
            if num_groups > max_groups_per_chrom:
                max_groups_per_chrom = num_groups
    
    print(f"所有染色体中最大的group数量: {max_groups_per_chrom}")
    
    # 计算每个染色体中最多的区间数（用于x轴范围）
    max_intervals_per_chrom = {}
    for chrom in chromosomes:
        if chrom in chrom_data:
            max_intervals = 0
            for group_data in chrom_data[chrom]['group_data'].values():
                if group_data['num_intervals'] > max_intervals:
                    max_intervals = group_data['num_intervals']
            max_intervals_per_chrom[chrom] = max_intervals
    
    # 统一的绘图参数
    bar_height = 0.8
    group_spacing = 0.4
    
    # 绘制每个染色体
    for idx, chrom in enumerate(chromosomes):
        if chrom not in chrom_data:
            continue
            
        data = chrom_data[chrom]
        target_groups = data['groups']
        group_data_dict = data['group_data']
        
        if not target_groups:
            continue
        
        # 计算子图位置
        row = idx // n_cols  # 0 或 1
        col = idx % n_cols
        
        # 调整行位置：第一行（chr1-chr6）正常，第二行（chr7-chr12）
        if row == 0:  # chr1-chr6
            ax = fig.add_subplot(gs[row * 2:row * 2 + 2, col])
        else:  # chr7-chr12
            # 这里的 index 3:5 紧跟在 Gap(2) 之后
            ax = fig.add_subplot(gs[row * 2 + 1:row * 2 + 3, col])
        
        # 获取该染色体中最多的区间数（用于归一化）
        max_intervals = max_intervals_per_chrom[chrom]
        
        # 绘制每个group
        for i, group in enumerate(target_groups):
            group_data = group_data_dict[group]
            haps = group_data['haps']
            num_intervals = group_data['num_intervals']
            
            y_pos = i * (bar_height + group_spacing)
            
            # 绘制每个区间
            for j in range(num_intervals):
                # 计算x位置（归一化到0-1，基于该染色体最大区间数）
                x_pos = j / max_intervals
                
                # 计算方块宽度
                width = 1 / max_intervals
                
                # 获取颜色（固定颜色，不使用score）
                color = get_color_for_hap(haps[j])
                
                # 绘制方块
                rect = plt.Rectangle((x_pos, y_pos), width, bar_height, 
                                   facecolor=color, edgecolor='black', linewidth=0.3)
                ax.add_patch(rect)
        
        # 设置统一的y轴范围
        y_max = max_groups_per_chrom * (bar_height + group_spacing)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, y_max)
        
        # 设置y轴标签（group名称）
        y_ticks = [i * (bar_height + group_spacing) + bar_height/2 for i in range(len(target_groups))]
        ax.set_yticks(y_ticks)
        ax.set_yticklabels(target_groups, fontsize=9, fontweight='bold')
        
        # 设置x轴
        labelpad = 20 if row == 1 else 15  # 修正判断条件
        ax.set_xlabel('Genomic Position (normalized)', fontsize=10, labelpad=labelpad)
        
        # 设置x轴刻度（基于该染色体最大区间数）
        xticks = np.linspace(0, 1, 6)
        # 计算实际长度（假设每个区间500kb）
        actual_length_mb = max_intervals * 0.5  # 转换为Mb
        xticklabels = [f'{i*actual_length_mb/5:.1f}Mb' for i in range(6)]
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticklabels, fontsize=8, rotation=45)
        
        # ================= 修改处 2 =================
        # 添加标题 - 调整chr7-chr12的标题位置
        # 之前的代码用 if col >= 6 判断是错的，因为col在0-5之间。
        # 这里改用 if row == 1 来判断是否为第二行。
        if row == 1:  # chr7-chr12
            title_y = 1.08  # 略微上调，使其紧凑但不过分远离
        else:  # chr1-chr6
            title_y = 1.05
        # ============================================
        
        ax.set_title(f'{chrom} ({len(target_groups)} groups)', 
                    fontsize=11, fontweight='bold', pad=15, y=title_y)
        
        # 添加背景网格
        ax.grid(True, alpha=0.2, linestyle=':', linewidth=0.5)
        
        # 对于group数量较少的染色体，添加空白标记
        if len(target_groups) < max_groups_per_chrom:
            # 在空白区域添加浅灰色背景
            empty_y_start = len(target_groups) * (bar_height + group_spacing)
            empty_height = y_max - empty_y_start
            if empty_height > 0:
                rect = plt.Rectangle((0, empty_y_start), 1, empty_height, 
                                   facecolor=(0.95, 0.95, 0.95), edgecolor='none', 
                                   alpha=0.3, zorder=0)
                ax.add_patch(rect)
    
    # 创建独立的图例区域（右侧）- 调整位置适应新的布局
    legend_ax = fig.add_axes([0.92, 0.25, 0.06, 0.4])
    legend_ax.axis('off')
    
    # 创建颜色图例（使用实际颜色块）
    legend_y = 0.95
    legend_spacing = 0.08
    
    # 图例标题
    legend_ax.text(0.5, legend_y, "Color Legend", 
                 fontsize=12, fontweight='bold', 
                 ha='center', transform=legend_ax.transAxes)
    
    # 绘制颜色块和标签
    for i, (label, color) in enumerate(zip(hap_labels, hap_rgb)):
        y_pos = legend_y - (i + 1) * legend_spacing
        
        # 绘制颜色块
        color_rect = plt.Rectangle((0.1, y_pos - 0.03), 0.2, 0.05, 
                                 facecolor=color, edgecolor='black', 
                                 transform=legend_ax.transAxes)
        legend_ax.add_patch(color_rect)
        
        # 添加标签
        legend_ax.text(0.35, y_pos - 0.005, label, 
                     fontsize=10, va='center',
                     transform=legend_ax.transAxes)
    
    # 添加说明文本
    info_y = legend_y - 7 * legend_spacing
    legend_ax.text(0.5, info_y, 
                 f"Each rectangle represents\na 500kb interval",
                 fontsize=9, ha='center', style='italic',
                 transform=legend_ax.transAxes)
    
    legend_ax.text(0.5, info_y - 0.08, 
                 f"All groups have same\nvertical height",
                 fontsize=9, ha='center', style='italic',
                 transform=legend_ax.transAxes)
    
    # 添加总标题
    plt.suptitle('Genomic Haplotig Assignment Heatmap by Chromosome', 
                fontsize=18, fontweight='bold', y=0.98)
    
    # 添加副标题
    fig.text(0.5, 0.94, 'Color indicates haplotig type. All groups have equal vertical height for comparison', 
            fontsize=14, ha='center', style='italic')
    
    # 添加底部说明
    fig.text(0.5, 0.02, 'X-axis normalized within each chromosome. Each chromosome shows its dominant groups.', 
            fontsize=10, ha='center', style='italic')
    
    # 调整布局
    # 稍微减小 hspace 从 0.5 到 0.4，配合 height_ratios 使用
    plt.subplots_adjust(left=0.05, right=0.9, bottom=0.05, top=0.85, 
                        hspace=0.4, wspace=0.3)
    
    # 保存图形
    output_path = os.path.join(output_dir, 'all_chromosomes_combined.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"已保存组合热力图: {output_path}")
    
    # 保存统计信息
    stats_path = os.path.join(output_dir, 'all_chromosomes_stats.txt')
    with open(stats_path, 'w') as f:
        f.write("Chromosome Distribution Statistics\n")
        f.write("=" * 50 + "\n\n")
        
        total_groups = 0
        for chrom in chromosomes:
            if chrom in chrom_data:
                num_groups = len(chrom_data[chrom]['groups'])
                total_groups += num_groups
                f.write(f"{chrom}: {num_groups} groups\n")
                f.write(f"  Groups: {', '.join(chrom_data[chrom]['groups'])}\n")
                
                # 每个group的区间数
                interval_counts = []
                for group in chrom_data[chrom]['groups']:
                    num_intervals = chrom_data[chrom]['group_data'][group]['num_intervals']
                    interval_counts.append(num_intervals)
                    f.write(f"    {group}: {num_intervals} intervals ({num_intervals*0.5:.1f} Mb)\n")
                
                if interval_counts:
                    min_intervals = min(interval_counts)
                    max_intervals = max(interval_counts)
                    f.write(f"    Interval range: {min_intervals}-{max_intervals} ({min_intervals*0.5:.1f}-{max_intervals*0.5:.1f} Mb)\n")
                f.write("\n")
        
        f.write(f"\nSummary Statistics:\n")
        f.write(f"- Total groups analyzed: {total_groups}\n")
        f.write(f"- Interval size: 500kb\n")
        f.write(f"- Color coding: red=hap1, blue=hap2, green=hap3, purple=hap4, brown=other_chrom, yellow=unreliable\n")
        f.write(f"- Note: All groups have equal vertical height for comparison\n")
        f.write(f"- X-axis normalized within each chromosome based on its longest group\n")
    
    print(f"已保存统计信息: {stats_path}")
    
    # 打印简要统计
    print(f"\n简要统计:")
    print(f"总group数: {total_groups}")
    for chrom in chromosomes:
        if chrom in chrom_data:
            groups = chrom_data[chrom]['groups']
            print(f"{chrom}: {len(groups)} groups")
            for group in groups:
                num_intervals = chrom_data[chrom]['group_data'][group]['num_intervals']
                print(f"  {group}: {num_intervals} intervals ({num_intervals*0.5:.1f} Mb)")

def main():
    # 设置数据目录
    data_dir = '/home/chenzh/HapHiC/data/C88/kmer/kmer_our'
    
    print("开始处理group文件...")
    
    # 处理所有group文件
    groups_data, group_dominant_chrom = process_all_groups(data_dir)
    
    print(f"\n处理完成！共找到 {len(groups_data)} 个groups")
    
    # 统计每个染色体的groups数量
    chrom_counts = defaultdict(int)
    for group, chrom in group_dominant_chrom.items():
        chrom_counts[chrom] += 1
    
    print("\n染色体分布统计:")
    for chrom in sorted(chrom_counts.keys()):
        print(f"  {chrom}: {chrom_counts[chrom]} groups")
    
    # 输出Bash格式的染色体分组信息
    print_chromosome_groups_bash_format(group_dominant_chrom)
    
    # 绘制组合热力图
    print("\n开始绘制组合热力图...")
    plot_all_chromosomes_combined(groups_data, group_dominant_chrom)
    
    print("\n组合热力图已生成！")

if __name__ == "__main__":
    # 设置matplotlib
    try:
        plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Helvetica', 'sans-serif']
        plt.rcParams['axes.unicode_minus'] = False
    except:
        pass
    
    main()