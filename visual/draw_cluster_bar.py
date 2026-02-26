import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from math import pi
import matplotlib.patches as mpatches

# 首先读取之前保存的CSV文件
csv_filename = 'clustering_results_updated.csv'
df = pd.read_csv(csv_filename)

# 定义指标列
metric_columns = ['Contiguity', 'Inter_homo_error_rate', 'Inter_nonhomo_error_rate', 
                 'total_error_rate', 'Anchoring_rate', 'Adjusted_contiguity']

# 指标显示名称（带单位）
metric_names = {
    'Contiguity': 'Contiguity',
    'Inter_homo_error_rate': 'Inter-homo Error Rate (%)',
    'Inter_nonhomo_error_rate': 'Inter-nonhomo Error Rate (%)',
    'total_error_rate': 'Total Error Rate (%)',
    'Anchoring_rate': 'Anchoring Rate (%)',
    'Adjusted_contiguity': 'Adjusted Contiguity (%)'
}

# 定义Paired配色方案（全局变量）
PAIRED_COLORS = ['#A6CEE3', '#1F78B4', '#B2DF8A', '#33A02C', '#FB9A99', '#E31A1C']
METHOD_COLORS = {
    'ConHiC': PAIRED_COLORS[0],  # 浅蓝
    'HapHiC': PAIRED_COLORS[1],    # 深蓝
    'ALLHiC': PAIRED_COLORS[2]     # 浅绿
}

# ========== 关键修改点：调整y轴起始百分比 ==========
# 设置Contiguity、Anchoring_rate、Adjusted_contiguity的y轴起始百分比
# 数值越大，差异显示越明显（90%意味着只显示顶部10%的区域）
Y_AXIS_START_PERCENTAGE = 0.90  # 从最小值的90%开始，让差异更明显
# ==============================================

def format_value(value, metric):
    """
    格式化数值，保留两位小数，避免出现0.0
    """
    if pd.isna(value):
        return '-'
    
    if metric == 'Contiguity':
        # Contiguity保留两位小数，但如果小于0.01则显示为<0.01
        if value < 0.01:
            return '<0.01'
        else:
            return f'{value:.2f}'
    else:
        # 其他指标保留两位小数，但如果小于0.01则显示为<0.01
        if value < 0.01:
            return '<0.01'
        else:
            return f'{value:.2f}'

def calculate_ylim_with_min_percentage(values, metric_type='general'):
    """
    计算y轴范围
    values: 数据值列表
    metric_type: 指标类型 
        - 'contiguity': Contiguity (从最小值的指定百分比开始)
        - 'percentage': 百分比指标 (从最小值的指定百分比开始)
        - 'error_rate': 错误率指标 (从0开始或从最小值的50%开始)
        - 'general': 通用情况
    """
    # 过滤掉NaN值
    valid_values = [v for v in values if not np.isnan(v)]
    if not valid_values:
        return 0, 1  # 默认范围
    
    min_val = min(valid_values)
    max_val = max(valid_values)
    
    # 根据指标类型设置y轴下限
    if metric_type == 'contiguity':
        # Contiguity在0-1之间，从最小值的指定百分比开始
        y_min = max(0, min_val * Y_AXIS_START_PERCENTAGE)  # 从最小值的指定百分比开始，但不小于0
        y_max = min(1.0, max_val * 1.15)  # 最大值向上留15%空间
        
        # 如果最小值非常接近0，则从0开始
        if min_val < 0.05:
            y_min = 0
            
    elif metric_type == 'percentage':
        # 百分比指标在0-100之间，从最小值的指定百分比开始
        y_min = max(0, min_val * Y_AXIS_START_PERCENTAGE)  # 从最小值的指定百分比开始，但不小于0
        y_max = min(100, max_val * 1.15)  # 最大值向上留15%空间
        
        # 如果最小值非常接近0，则从0开始
        if min_val < 5:
            y_min = 0
            
    elif metric_type == 'error_rate':
        # 错误率指标
        if min_val < 1:  # 如果最小值很小，从0开始更合适
            y_min = 0
        else:
            y_min = min_val * 0.5  # 从最小值的一半开始
        y_max = max_val * 1.2  # 最大值向上留20%空间
    else:
        # 通用情况
        y_min = min_val * 0.5
        y_max = max_val * 1.2
    
    return y_min, y_max

# 创建综合数据展示的柱状图（保持原文件名）
def create_comprehensive_bar_chart(df, metric_columns, metric_names):
    """创建包含6个子图的综合柱状图，每个子图对应一个指标（保持原文件名）"""
    
    # 准备数据：将C88的2M和0.8M分开处理
    data_groups = []
    group_labels = []
    
    # 处理C88
    c88_2m = df[(df['Species'] == 'C88') & (df['N50'] == 'N50=2M')]
    if not c88_2m.empty:
        data_groups.append(c88_2m)
        group_labels.append('C88 (2M)')
    
    c88_08m = df[(df['Species'] == 'C88') & (df['N50'] == 'N50=0.8M')]
    if not c88_08m.empty:
        data_groups.append(c88_08m)
        group_labels.append('C88 (0.8M)')
    
    # 处理其他物种
    other_species = [s for s in df['Species'].unique() if s != 'C88']
    for species in other_species:
        species_df = df[df['Species'] == species]
        data_groups.append(species_df)
        group_labels.append(species)
    
    # 定义方法列表
    methods = ['ConHiC', 'HapHiC', 'ALLHiC']
    
    # 创建图形 - 2行3列，每个子图更大
    fig, axes = plt.subplots(2, 3, figsize=(24, 15))
    axes = axes.flatten()
    
    # 为每个指标创建子图
    for idx, metric in enumerate(metric_columns):
        ax = axes[idx]
        
        # 准备x轴位置
        n_groups = len(data_groups)
        x = np.arange(n_groups)  # 每个组的位置
        width = 0.25  # 柱子的宽度
        
        # 存储每个方法在每个组中的数值
        method_values = {method: [] for method in methods}
        
        # 收集数据
        for group_df in data_groups:
            for method in methods:
                method_data = group_df[group_df['Method'] == method]
                if not method_data.empty and metric in method_data.columns:
                    value = method_data[metric].iloc[0]
                    if not pd.isna(value):
                        method_values[method].append(value)
                    else:
                        method_values[method].append(np.nan)
                else:
                    method_values[method].append(np.nan)
        
        # 绘制柱状图
        bars = []
        for i, method in enumerate(methods):
            # 计算每个方法柱子的x位置
            x_pos = x + (i - 1) * width  # 居中排列：-1, 0, 1 * width
            
            # 绘制柱子
            bar = ax.bar(x_pos, method_values[method], width, 
                        label=method, color=METHOD_COLORS[method], alpha=0.8,
                        edgecolor='black', linewidth=0.5)
            bars.append(bar)
        
        # 收集所有有效值用于计算y轴范围
        all_values = []
        for method in methods:
            all_values.extend([v for v in method_values[method] if not np.isnan(v)])
        
        # 根据指标类型设置y轴范围
        if metric == 'Contiguity':
            y_min, y_max = calculate_ylim_with_min_percentage(all_values, 'contiguity')
        elif metric in ['Anchoring_rate', 'Adjusted_contiguity']:
            y_min, y_max = calculate_ylim_with_min_percentage(all_values, 'percentage')
        elif metric in ['total_error_rate', 'Inter_homo_error_rate', 'Inter_nonhomo_error_rate']:
            y_min, y_max = calculate_ylim_with_min_percentage(all_values, 'error_rate')
        else:
            y_min, y_max = calculate_ylim_with_min_percentage(all_values, 'general')
        
        ax.set_ylim(y_min, y_max)
        
        # 添加数值标签（保留两位小数）
        for i, method in enumerate(methods):
            x_pos = x + (i - 1) * width
            values = method_values[method]
            
            for j, (pos, val) in enumerate(zip(x_pos, values)):
                if not np.isnan(val):
                    # 格式化标签
                    label = format_value(val, metric)
                    
                    # 添加标签（在柱子顶部）
                    ax.text(pos, val + (y_max - y_min) * 0.02, 
                           label, ha='center', va='bottom', fontsize=8, rotation=45)
        
        # 设置x轴标签
        ax.set_xticks(x)
        ax.set_xticklabels(group_labels, rotation=45, ha='right', fontsize=10)
        
        # 设置y轴标签
        if metric in ['Contiguity']:
            ax.set_ylabel('Value', fontsize=11)
        else:
            ax.set_ylabel('Percentage (%)', fontsize=11)
        
        # 设置标题 - 只保留子图标题
        ax.set_title(metric_names[metric], fontsize=13, fontweight='bold', pad=10)
        
        # 添加网格线
        ax.grid(True, alpha=0.3, axis='y', linestyle='--')
        
        # 在y轴上标记最小值的位置（可选，用于可视化参考）
        if len(all_values) > 0:
            min_val = min(all_values)
            if metric in ['Contiguity', 'Anchoring_rate', 'Adjusted_contiguity']:
                ax.axhline(y=min_val * Y_AXIS_START_PERCENTAGE, color='gray', linestyle=':', alpha=0.5, linewidth=0.8, 
                          label=f'{int(Y_AXIS_START_PERCENTAGE*100)}% of min')
        
        # 移除顶部和右侧边框
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        # 加粗左侧和底部边框
        ax.spines['left'].set_linewidth(1.5)
        ax.spines['bottom'].set_linewidth(1.5)
    
    # 创建图例 - 放在图表顶部
    legend_elements = [mpatches.Patch(color=METHOD_COLORS[method], label=method, alpha=0.8, edgecolor='black') 
                      for method in methods]
    
    # 调整图例位置
    fig.legend(handles=legend_elements, loc='upper center', 
               bbox_to_anchor=(0.5, 1.02), ncol=3, fontsize=13, 
               frameon=True, framealpha=0.9, edgecolor='gray')
    
    # 完全移除了总标题
    
    # 调整布局，给图例留出空间
    plt.tight_layout(rect=[0, 0, 1, 0.97])
    
    # 保存图像 - 保持原文件名不变
    output_filename = 'clustering_results_table_compact.png'
    plt.savefig(output_filename, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    
    print(f"综合柱状图已保存为: {output_filename} (已移除顶部大标题，数值保留两位小数，Contiguity/Anchoring/Adjusted从最小值的{int(Y_AXIS_START_PERCENTAGE*100)}%开始)")
    
    return output_filename

# 创建数值汇总表格
def create_value_summary_table(df, metric_columns):
    """创建数值汇总表格，显示所有数值（保留两位小数）"""
    
    print("\n" + "="*80)
    print("数值汇总表格（所有数值保留两位小数）:")
    
    # 准备数据
    summary_data = []
    
    # 处理C88
    c88_2m = df[(df['Species'] == 'C88') & (df['N50'] == 'N50=2M')]
    if not c88_2m.empty:
        for _, row in c88_2m.iterrows():
            summary_data.append(('C88 (2M)', row['Method'], row))
    
    c88_08m = df[(df['Species'] == 'C88') & (df['N50'] == 'N50=0.8M')]
    if not c88_08m.empty:
        for _, row in c88_08m.iterrows():
            summary_data.append(('C88 (0.8M)', row['Method'], row))
    
    # 处理其他物种
    other_species = [s for s in df['Species'].unique() if s != 'C88']
    for species in other_species:
        species_df = df[df['Species'] == species]
        for _, row in species_df.iterrows():
            summary_data.append((species, row['Method'], row))
    
    # 打印表格
    print("\n{:<15} {:<10} {:<12} {:<15} {:<18} {:<15} {:<15} {:<18}".format(
        "Species", "Method", "Contiguity", "HomoErr%", "NonHomoErr%", "TotalErr%", "Anchor%", "AdjContig%"))
    print("-" * 130)
    
    for species, method, row in summary_data:
        cont = format_value(row['Contiguity'], 'Contiguity')
        homo_err = format_value(row['Inter_homo_error_rate'], 'error')
        nonhomo_err = format_value(row['Inter_nonhomo_error_rate'], 'error')
        total_err = format_value(row['total_error_rate'], 'error')
        anchor = format_value(row['Anchoring_rate'], 'percentage')
        adj_cont = format_value(row['Adjusted_contiguity'], 'percentage')
        
        print("{:<15} {:<10} {:<12} {:<15} {:<18} {:<15} {:<15} {:<18}".format(
            species, method, cont, homo_err, nonhomo_err, total_err, anchor, adj_cont))
    
    # 保存到CSV
    summary_rows = []
    for species, method, row in summary_data:
        summary_rows.append({
            'Species': species,
            'Method': method,
            'Contiguity': format_value(row['Contiguity'], 'Contiguity'),
            'Inter_homo_error_rate': format_value(row['Inter_homo_error_rate'], 'error'),
            'Inter_nonhomo_error_rate': format_value(row['Inter_nonhomo_error_rate'], 'error'),
            'total_error_rate': format_value(row['total_error_rate'], 'error'),
            'Anchoring_rate': format_value(row['Anchoring_rate'], 'percentage'),
            'Adjusted_contiguity': format_value(row['Adjusted_contiguity'], 'percentage')
        })
    
    summary_df = pd.DataFrame(summary_rows)
    summary_csv = 'value_summary_2decimals.csv'
    summary_df.to_csv(summary_csv, index=False, encoding='utf-8-sig')
    print(f"\n数值汇总表格已保存为: {summary_csv}")

# 执行主函数
print("="*80)
print("开始创建综合柱状图（使用Paired配色方案，已移除顶部大标题）...")
print(f"配色方案:")
print(f"  - ConHiC: {METHOD_COLORS['ConHiC']} (浅蓝)")
print(f"  - HapHiC: {METHOD_COLORS['HapHiC']} (深蓝)")
print(f"  - ALLHiC: {METHOD_COLORS['ALLHiC']} (浅绿)")
print(f"\n数值格式: 所有数值保留两位小数，小于0.01的显示为'<0.01'")
print(f"y轴设置: Contiguity/Anchoring/Adjusted从最小值的{int(Y_AXIS_START_PERCENTAGE*100)}%开始（差异显示更明显）")
print(f"         错误率指标从0或最小值的一半开始")

# 只创建综合柱状图（不创建各物种单独柱状图）
comprehensive_chart = create_comprehensive_bar_chart(df, metric_columns, metric_names)

# 创建数值汇总表格
create_value_summary_table(df, metric_columns)

# 创建方法对比的汇总表（显示每个指标的最佳方法）
print("\n" + "="*80)
print("各指标最佳表现统计（数值保留两位小数）:")

# 定义数据组
data_groups = []
group_labels = []

# 处理C88
c88_2m = df[(df['Species'] == 'C88') & (df['N50'] == 'N50=2M')]
if not c88_2m.empty:
    data_groups.append(('C88 (2M)', c88_2m))

c88_08m = df[(df['Species'] == 'C88') & (df['N50'] == 'N50=0.8M')]
if not c88_08m.empty:
    data_groups.append(('C88 (0.8M)', c88_08m))

# 处理其他物种
other_species = [s for s in df['Species'].unique() if s != 'C88']
for species in other_species:
    species_df = df[df['Species'] == species]
    data_groups.append((species, species_df))

# 为每个指标找出最佳方法
for metric in metric_columns:
    print(f"\n{metric_names[metric]}:")
    
    for group_name, group_df in data_groups:
        # 过滤NaN值
        valid_data = group_df[~group_df[metric].isna()].copy()
        if len(valid_data) == 0:
            continue
        
        # 确定最佳值方向
        if metric in ['total_error_rate', 'Inter_homo_error_rate', 'Inter_nonhomo_error_rate']:
            # 错误率越低越好
            best_row = valid_data.loc[valid_data[metric].idxmin()]
            best_value = best_row[metric]
            best_method = best_row['Method']
            symbol = "↓"  # 越低越好
        else:
            # 其他指标越高越好
            best_row = valid_data.loc[valid_data[metric].idxmax()]
            best_value = best_row[metric]
            best_method = best_row['Method']
            symbol = "↑"  # 越高越好
        
        # 格式化输出
        formatted_value = format_value(best_value, metric)
        print(f"  {group_name}: {best_method} ({formatted_value}) {symbol}")

# 创建方法胜率统计
print("\n" + "="*80)
print("方法胜率统计（每个分组内获得最佳表现的次数）:")

# 统计每个方法在每组中获胜的次数
method_wins = {'ConHiC': 0, 'HapHiC': 0, 'ALLHiC': 0}
total_comparisons = 0

for group_name, group_df in data_groups:
    for metric in metric_columns:
        valid_data = group_df[~group_df[metric].isna()].copy()
        if len(valid_data) == 0:
            continue
        
        total_comparisons += 1
        
        if metric in ['total_error_rate', 'Inter_homo_error_rate', 'Inter_nonhomo_error_rate']:
            best_method = valid_data.loc[valid_data[metric].idxmin(), 'Method']
        else:
            best_method = valid_data.loc[valid_data[metric].idxmax(), 'Method']
        
        method_wins[best_method] += 1

# 计算胜率
print(f"\n总比较次数: {total_comparisons}")
for method in method_wins:
    win_rate = (method_wins[method] / total_comparisons) * 100
    print(f"  {method}: {method_wins[method]} 次获胜 ({win_rate:.1f}%)")

# 创建方法综合评分
print("\n" + "="*80)
print("方法综合评分（基于排名加权）:")

# 为每个方法计算综合得分
method_scores = {'ConHiC': 0, 'HapHiC': 0, 'ALLHiC': 0}
total_rankings = 0

for group_name, group_df in data_groups:
    for metric in metric_columns:
        valid_data = group_df[~group_df[metric].isna()].copy()
        if len(valid_data) == 0:
            continue
        
        total_rankings += 1
        
        # 确定排序方向
        if metric in ['total_error_rate', 'Inter_homo_error_rate', 'Inter_nonhomo_error_rate']:
            sorted_data = valid_data.sort_values(by=metric, ascending=True)
        else:
            sorted_data = valid_data.sort_values(by=metric, ascending=False)
        
        # 给排名赋分（第一名3分，第二名2分，第三名1分）
        for rank, (idx, row) in enumerate(sorted_data.iterrows()):
            method = row['Method']
            score = 3 - rank  # 第一名3分，第二名2分，第三名1分
            method_scores[method] += score

# 计算平均分
print("\n各方法平均得分（最高3分）:")
for method in method_scores:
    avg_score = method_scores[method] / total_rankings if total_rankings > 0 else 0
    print(f"  {method}: {avg_score:.2f} 分")

print("\n" + "="*80)
print("所有图表生成完成!")
print(f"主要输出文件:")
print(f"  1. 综合柱状图: {comprehensive_chart} (已移除顶部大标题，文件名保持不变)")
print(f"  2. 数值汇总表格: value_summary_2decimals.csv")
print(f"\n主要改进:")
print(f"  - 已移除综合图的顶部大标题")
print(f"  - 所有数值保留两位小数，小于0.01的显示为'<0.01'")
print(f"  - Contiguity、Anchoring rate、Adjusted contiguity的y轴从最小值的{int(Y_AXIS_START_PERCENTAGE*100)}%开始（差异更明显）")
print(f"  - 错误率指标根据数值大小智能调整y轴起始点")
print(f"  - 添加灰色虚线标记最小值{int(Y_AXIS_START_PERCENTAGE*100)}%的位置（参考线）")
print(f"  - 生成了数值汇总表格，方便查看所有数据")