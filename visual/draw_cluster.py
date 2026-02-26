import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import re
from math import pi

# 更新后的原始数据（包含CS物种的更正）
data = """C88						
N50=2M	Contiguity	Inter_homo_error_rate	Inter_nonhomo_error_rate	total_error_rate	Anchoring rate	Adjusted contiguity
ALLHiC	0.8651 	11.6552 	4.9323 	16.5875 	82.2149 	71.1307 
HapHiC	0.8472 	6.0354 	5.4297 	11.4651 	88.5367 	75.0082 
ConHiC	0.8732 	5.4975 	4.9570 	10.4545 	89.9603 	78.5582 
N50=0.8M	Contiguity	Inter_homo_error_rate	Inter_nonhomo_error_rate	total_error_rate	Anchoring rate	Adjusted contiguity
ALLHiC	0.8461 	14.6235 	3.2481 	17.8716 	77.8879 	65.9010 
HapHiC	0.8484 	6.7157 	5.5210 	12.2367 	83.5571 	70.8898 
ConHiC	0.8531 	6.6996 	5.4735 	12.1731 	85.9468 	73.3212 
XinJiangDaYe						
	Contiguity	Inter_homo_error_rate	Inter_nonhomo_error_rate	total_error_rate	Anchoring rate	Adjusted contiguity
ALLHiC	0.5213 	65.0584 	0.7657 	65.8241 	82.2385 	42.8709 
HapHiC	0.8695 	12.9909 	0.0462 	13.0371 	80.7125 	70.1795 
ConHiC	0.8706 	12.8872 	0.0462 	12.9334 	80.7117 	70.2676 
Saccharum_spontaneum_Np-X						
	Contiguity	Inter_homo_error_rate	Inter_nonhomo_error_rate	total_error_rate	Anchoring rate	Adjusted contiguity
ALLHiC	0.8267 	20.5660 	1.2184 	21.7844 	98.6683 	81.5691 
HapHiC	0.9096 	11.0376 	1.0788 	12.1164 	90.5987 	82.4086 
ConHiC	0.9323 	10.6879 	0.9990 	11.6869 	90.4177 	84.2964 
ZM-4						
	Contiguity	Inter_homo_error_rate	Inter_nonhomo_error_rate	total_error_rate	Anchoring rate	Adjusted contiguity
ALLHiC	0.6793 	32.4232 	7.9408 	40.3640 	92.5009 	62.8359 
HapHiC	0.8221 	12.0228 	5.9876 	18.0104 	89.4927 	73.5719 
ConHiC	0.8369 	11.4906 	5.5193 	17.0099 	89.6777 	75.0513 
CS						
	Contiguity	Inter_homo_error_rate	Inter_nonhomo_error_rate	total_error_rate	Anchoring rate	Adjusted contiguity
ALLHiC	0.8788 	42.4796 	7.7443 	50.2239 	80.5699 	70.8048 
HapHiC	0.5585 	37.4164 	6.1224 	43.5388 	78.5898 	43.8924 
ConHiC	0.5587 	37.3903 	6.1188 	43.5091 	78.6272 	43.9290 
AP						
	Contiguity	Inter_homo_error_rate	Inter_nonhomo_error_rate	total_error_rate	Anchoring rate	Adjusted contiguity
HapHiC	0.8288 	13.0630 	4.5476 	17.6106 	89.5684 	74.2343 
ConHiC	0.8276 	13.1345 	4.6200 	17.7545 	90.0634 	74.5365"""

# 简化的解析逻辑
lines = data.strip().split('\n')
records = []
current_species = None
current_n50 = "Default"
current_columns = None

i = 0
while i < len(lines):
    line = lines[i].strip()
    
    # 跳过空行
    if not line:
        i += 1
        continue
    
    # 检查是否为物种标题（不以制表符开头且不包含等号和常见指标名）
    parts = line.split('\t')
    if (len(parts) == 1 and 
        'N50=' not in line and 
        'Contiguity' not in line and 
        not any(method in line for method in ['ALLHiC', 'HapHiC', 'ConHiC'])):
        current_species = line
        current_n50 = "Default"
        i += 1
        continue
    
    # 检查是否为N50行
    if 'N50=' in line:
        current_n50 = parts[0].strip()
        current_columns = ['Method'] + [col.strip() for col in parts[1:] if col.strip()]
        i += 1
        continue
    
    # 检查是否为表头行（包含Contiguity）
    if 'Contiguity' in line:
        parts = [p.strip() for p in line.split('\t') if p.strip()]
        if len(parts) >= 6:  # 完整的表头
            current_columns = ['Method'] + parts
        i += 1
        continue
    
    # 检查是否为数据行（以方法名开头）
    if any(line.startswith(method) or f"\t{method}\t" in line for method in ['ALLHiC', 'HapHiC', 'ConHiC']):
        parts = [p.strip() for p in line.split('\t') if p.strip()]
        
        # 处理数据缺失的情况
        if len(parts) < 7 and len(parts) > 1:
            # 尝试处理不完整的数据行
            method = parts[0]
            values = parts[1:]
            
            # 创建记录
            record = {
                'Species': current_species,
                'N50': current_n50,
                'Method': method
            }
            
            # 使用固定的列名添加指标值
            fixed_cols = ['Contiguity', 'Inter_homo_error_rate', 'Inter_nonhomo_error_rate', 
                         'total_error_rate', 'Anchoring_rate', 'Adjusted_contiguity']
            
            for idx in range(6):  # 尝试填充6个指标
                if idx < len(values) and values[idx].strip():
                    try:
                        record[fixed_cols[idx]] = float(values[idx])
                    except ValueError:
                        record[fixed_cols[idx]] = np.nan
                else:
                    record[fixed_cols[idx]] = np.nan
            
            records.append(record)
        
        elif len(parts) >= 7:  # 完整的数据行
            method = parts[0]
            values = parts[1:]
            
            # 创建记录
            record = {
                'Species': current_species,
                'N50': current_n50,
                'Method': method
            }
            
            # 使用固定的列名添加指标值
            fixed_cols = ['Contiguity', 'Inter_homo_error_rate', 'Inter_nonhomo_error_rate', 
                         'total_error_rate', 'Anchoring_rate', 'Adjusted_contiguity']
            
            for idx, (col, val) in enumerate(zip(fixed_cols, values[:6])):  # 只取前6个值
                try:
                    record[col] = float(val)
                except ValueError:
                    record[col] = np.nan  # 对于无法转换的值，使用NaN
            
            records.append(record)
    
    i += 1

# 创建DataFrame
df = pd.DataFrame(records)

# 简化和重命名物种名称
species_mapping = {
    'C88': 'C88',
    'XinJiangDaYe': 'XJDY',
    'Saccharum_spontaneum_Np-X': 'SS',
    'ZM-4': 'ZM-4',
    'CS': 'CS',
    'AP': 'AP'
}
df['Species'] = df['Species'].map(species_mapping)

# 重新排列列顺序
required_columns = ['Species', 'N50', 'Method']
other_columns = [col for col in df.columns if col not in required_columns]
column_order = required_columns + other_columns
df = df[column_order]

# 显示数据概览
print("="*80)
print("数据概览:")
print(f"总行数: {len(df)}")
print(f"物种列表: {df['Species'].unique().tolist()}")
print(f"方法列表: {sorted(df['Method'].unique().tolist())}")
print(f"N50列表: {sorted(df['N50'].unique().tolist())}")

# 保存为CSV文件
csv_filename = 'clustering_results_updated.csv'
df.to_csv(csv_filename, index=False, encoding='utf-8-sig')
print(f"\n数据已保存为CSV文件: {csv_filename}")

# 显示整理后的数据
print("\n整理后的数据表格:")
print(df.to_string())

# 定义分组函数
def get_groups_with_spacing(df):
    """将数据分成不同的组，C88的2M和0.8M分开处理，并添加分组信息"""
    groups = []
    
    # C88的特殊处理：2M和0.8M分开
    if 'C88' in df['Species'].values:
        # C88 N50=2M
        c88_2m = df[(df['Species'] == 'C88') & (df['N50'] == 'N50=2M')]
        if not c88_2m.empty:
            groups.append(('C88 (N50=2M)', c88_2m, True))
            
        # C88 N50=0.8M
        c88_08m = df[(df['Species'] == 'C88') & (df['N50'] == 'N50=0.8M')]
        if not c88_08m.empty:
            groups.append(('C88 (N50=0.8M)', c88_08m, True))
        
        # 其他物种
        other_species = df[~((df['Species'] == 'C88') & (df['N50'].isin(['N50=2M', 'N50=0.8M'])))]
        for species in sorted(other_species['Species'].unique()):
            species_df = other_species[other_species['Species'] == species]
            if not species_df.empty:
                groups.append((f'{species}', species_df, True))
    else:
        # 如果没有C88，按物种分组
        for species in sorted(df['Species'].unique()):
            species_df = df[df['Species'] == species]
            groups.append((species, species_df, True))
    
    return groups

# 获取分组
groups = get_groups_with_spacing(df)

# 创建美观的表格图像，并将最佳值加粗，添加分组间隔
def create_table_image_with_best_values_and_spacing(df, groups, filename='clustering_results_table_compact.png'):
    # 设置中文字体（如果系统中存在）
    try:
        plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans']
        plt.rcParams['axes.unicode_minus'] = False
    except:
        pass
    
    # 调整图形尺寸，缩小列宽
    fig, ax = plt.subplots(figsize=(22, 18))  # 增加宽度，减少列间空隙
    ax.axis('tight')
    ax.axis('off')
    
    # 准备表格数据 - 包含间隔行
    table_data = []
    
    # 添加表头
    headers = df.columns.tolist()
    
    # 简化列名显示 - 更短的缩写
    display_headers = []
    column_mapping = {}
    column_display_names = {
        'Species': 'Species',
        'N50': 'N50',
        'Method': 'Method',
        'Contiguity': 'Contig',
        'Inter_homo_error_rate': 'HomoErr%',
        'Inter_nonhomo_error_rate': 'NonHomoErr%',
        'total_error_rate': 'TotalErr%',
        'Anchoring_rate': 'Anchor%',
        'Adjusted_contiguity': 'AdjContig%'
    }
    
    for h in headers:
        display_name = column_display_names.get(h, h)
        display_headers.append(display_name)
        column_mapping[display_name] = h
    
    table_data.append(display_headers)
    
    # 创建一个映射，将原始df行索引映射到表格行号
    df_to_table_row = {}
    table_row_counter = 1  # 从1开始（0是表头）
    
    # 标识每个单元格的格式（是否加粗）
    cell_formats = []
    # 初始化格式矩阵，所有单元格默认不加粗
    cell_formats.append([{'weight': 'normal'} for _ in range(len(headers))])
    
    # 定义要评估的指标列
    metric_columns = ['Contiguity', 'Inter_homo_error_rate', 'Inter_nonhomo_error_rate', 
                     'total_error_rate', 'Anchoring_rate', 'Adjusted_contiguity']
    
    # 存储最佳值信息用于后续标记
    best_values_info = []
    
    # 第一步：找出所有分组的最佳值（忽略NaN值）
    for group_idx, (group_name, group_df, is_new_group) in enumerate(groups):
        group_best_values = {}
        
        # 对于每个指标，找出最佳值
        for metric in metric_columns:
            if metric in group_df.columns:
                # 过滤掉NaN值
                valid_data = group_df[~group_df[metric].isna()]
                if len(valid_data) == 0:
                    continue
                
                # 确定优化方向
                if metric in ['total_error_rate', 'Inter_homo_error_rate', 'Inter_nonhomo_error_rate']:
                    # 错误率：越低越好
                    best_value = valid_data[metric].min()
                    best_indices = valid_data[valid_data[metric] == best_value].index.tolist()
                else:
                    # 其他指标：越高越好
                    best_value = valid_data[metric].max()
                    best_indices = valid_data[valid_data[metric] == best_value].index.tolist()
                
                # 对于并列第一的情况，只选择第一个
                best_index = best_indices[0]
                group_best_values[metric] = (best_index, best_value)
        
        best_values_info.append((group_name, group_best_values))
    
    # 第二步：添加数据行和间隔行，并标记最佳值
    first_group = True
    for group_idx, (group_name, group_df, is_new_group) in enumerate(groups):
        # 如果不是第一个分组，添加一个空行作为间隔
        if not first_group:
            # 添加间隔行
            table_data.append([''] * len(headers))
            cell_formats.append([{'weight': 'normal'} for _ in range(len(headers))])
            table_row_counter += 1
        
        first_group = False
        
        # 添加分组的所有数据行
        for idx, row in group_df.iterrows():
            # 格式化数值，统一保留4位小数，NaN显示为'-'
            formatted_row = []
            for val in row:
                if isinstance(val, (int, float, np.number)):
                    if pd.isna(val):
                        formatted_row.append('-')
                    else:
                        # 根据列名决定显示精度
                        col_name = headers[len(formatted_row)]
                        if col_name in ['Contiguity']:
                            formatted_row.append(f"{val:.4f}")
                        elif col_name in ['Inter_homo_error_rate', 'Inter_nonhomo_error_rate', 
                                         'total_error_rate', 'Anchoring_rate', 'Adjusted_contiguity']:
                            formatted_row.append(f"{val:.2f}")
                        else:
                            formatted_row.append(f"{val:.4f}")
                else:
                    formatted_row.append(str(val))
            table_data.append(formatted_row)
            
            # 记录映射关系
            df_to_table_row[idx] = table_row_counter
            table_row_counter += 1
            
            # 添加格式行
            cell_formats.append([{'weight': 'normal'} for _ in range(len(headers))])
    
    # 第三步：标记最佳值（只加粗，不标红）
    for group_idx, (group_name, group_best_values) in enumerate(best_values_info):
        for metric, (best_idx, best_value) in group_best_values.items():
            if best_idx in df_to_table_row:
                table_row = df_to_table_row[best_idx]
                
                # 找到这个指标在表格中的列号
                display_metric = None
                for disp, orig in column_mapping.items():
                    if orig == metric:
                        display_metric = disp
                        break
                
                if display_metric and display_metric in display_headers:
                    table_col = display_headers.index(display_metric)
                    # 只加粗，不标红
                    cell_formats[table_row][table_col]['weight'] = 'bold'
    
    # 创建表格
    table = ax.table(cellText=table_data, loc='center', cellLoc='center')
    
    # 设置表格样式 - 增大字号，调整缩放
    table.auto_set_font_size(False)
    table.set_fontsize(10)  # 增大到10号字
    table.scale(1.0, 1)   # 减少水平缩放，增加垂直缩放
    
    # 设置表头样式
    for i in range(len(headers)):
        table[(0, i)].set_facecolor('#2E75B6')  # 深蓝色
        table[(0, i)].set_text_props(weight='bold', color='white', fontsize=11)  # 增大表头字体
        table[(0, i)].set_height(0.025)  # 表头高度
    
    # 设置行颜色并应用格式
    for i in range(1, len(table_data)):
        # 检查是否是间隔行（所有单元格都为空）
        is_spacing_row = all(cell == '' for cell in table_data[i])
        
        if is_spacing_row:
            # 间隔行使用浅灰色背景
            for j in range(len(headers)):
                table[(i, j)].set_facecolor('#F5F5F5')
                table[(i, j)].set_height(0.015)  # 间隔行高度
        else:
            # 数据行交替颜色
            if i % 2 == 1:
                row_color = '#FFFFFF'  # 白色
            else:
                row_color = '#F9F9F9'  # 浅灰色
            
            for j in range(len(headers)):
                table[(i, j)].set_facecolor(row_color)
                
                # 应用格式（加粗）
                fmt = cell_formats[i][j]
                weight = fmt.get('weight', 'normal')
                table[(i, j)].set_text_props(weight=weight, fontsize=10)
                
                # 设置数据行高度
                table[(i, j)].set_height(0.02)
    
    # 设置Species、N50和Method列的特殊样式
    for i in range(1, len(table_data)):
        # 跳过间隔行
        if all(cell == '' for cell in table_data[i]):
            continue
            
        table[(i, 0)].set_facecolor('#E6F2FF')  # Species列 - 淡蓝色
        table[(i, 1)].set_facecolor('#FFF2E6')  # N50列 - 淡橙色
        table[(i, 2)].set_facecolor('#F0FFE6')  # Method列 - 淡绿色
    
    # 设置单元格边框
    for i in range(len(table_data)):
        for j in range(len(headers)):
            table[(i, j)].set_edgecolor('#DDDDDD')  # 浅灰色边框
            table[(i, j)].set_linewidth(0.5)  # 细边框
    
    # 自动调整列宽 - 根据内容调整
    # 计算每列的最大字符数
    col_widths = []
    for j in range(len(headers)):
        max_len = 0
        for i in range(len(table_data)):
            cell_text = str(table_data[i][j])
            max_len = max(max_len, len(cell_text))
        # 根据最大字符数设置列宽比例
        col_widths.append(max_len * 0.015)  # 调整系数以控制列宽
    
    # 应用列宽调整
    for j in range(len(headers)):
        table.auto_set_column_width([j])
    
    # 关键修改：将表格定位到图形边缘，完全去除留白
    ax.set_position([0, 0, 1, 1])  # 表格占据整个图形区域
    
    # 移除所有边距
    plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
    
    # 保存图像，设置bbox_inches为'tight'但不留白
    plt.savefig(filename, dpi=300, bbox_inches='tight', pad_inches=0, facecolor='white')
    print(f"\n紧凑表格图像已保存为: {filename}")
    
    return best_values_info

# 创建表格图像
best_values_info = create_table_image_with_best_values_and_spacing(df, groups)

# ============================================================================
# 改进的条形图绘制函数
# ============================================================================

def create_improved_bar_charts(df):
    """创建改进的标准化条形图，使用原始数值，优化显示，y轴从最小数值的一半开始"""
    
    # 准备数据
    species_list = []
    for species in df['Species'].unique():
        if species == 'C88':
            for n50_val in ['N50=2M', 'N50=0.8M']:
                species_list.append(f"{species} ({n50_val})")
        else:
            species_list.append(species)
    
    metrics = ['Contiguity', 'Anchoring_rate', 'Adjusted_contiguity',
               'Inter_homo_error_rate', 'Inter_nonhomo_error_rate', 'total_error_rate']
    
    # 指标显示名称（带单位）
    metric_names = {
        'Contiguity': 'Contiguity',
        'Anchoring_rate': 'Anchoring Rate (%)',
        'Adjusted_contiguity': 'Adjusted Contiguity (%)',
        'Inter_homo_error_rate': 'Homo Error Rate (%)',
        'Inter_nonhomo_error_rate': 'Non-homo Error Rate (%)',
        'total_error_rate': 'Total Error Rate (%)'
    }
    
    # 为每个物种创建条形图
    for species_full in species_list:
        # 解析物种名称
        if "C88" in species_full:
            if "(N50=2M)" in species_full:
                species = 'C88'
                n50_val = 'N50=2M'
                species_df = df[(df['Species'] == species) & (df['N50'] == n50_val)]
            else:
                species = 'C88'
                n50_val = 'N50=0.8M'
                species_df = df[(df['Species'] == species) & (df['N50'] == n50_val)]
        else:
            species = species_full
            species_df = df[df['Species'] == species]
        
        if species_df.empty:
            continue
        
        # 创建图形 - 调整为2行3列，每个子图更大
        fig, axes = plt.subplots(2, 3, figsize=(20, 14))
        axes = axes.flatten()
        
        # 颜色方案 - 统一使用HapHiC
        colors = {
            'ALLHiC': '#4ECDC4',   # 青色
            'HapHiC': '#45B7D1',   # 蓝色
            'ConHiC': '#FF6B6B'   # 红色
        }
        
        # 为每个指标创建子图
        for idx, metric in enumerate(metrics):
            ax = axes[idx]
            
            # 获取数据（过滤掉NaN值）
            methods = []
            values = []
            for method in species_df['Method'].unique():
                method_data = species_df[species_df['Method'] == method]
                if metric in method_data.columns and not pd.isna(method_data[metric].iloc[0]):
                    methods.append(method)
                    value = method_data[metric].iloc[0]
                    values.append(value)
            
            # 如果没有有效数据，跳过这个指标
            if len(values) == 0:
                ax.axis('off')
                continue
            
            # 计算y轴范围 - 从最小值的一半开始
            max_value = max(values)
            min_value = min(values)
            
            # 确定y轴范围
            if metric in ['Contiguity', 'Anchoring_rate', 'Adjusted_contiguity']:
                # 对于0-100范围的指标
                if metric == 'Contiguity':
                    # Contiguity在0-1之间
                    y_min = max(0, min_value * 0.5)  # 从最小值的一半开始，但不小于0
                    y_max = min(1.0, max_value * 1.25)  # 最大值向上留25%空间
                    if max_value < 0.9:
                        y_max = max_value * 1.35  # 如果最大值较小，留更多空间
                else:
                    # Anchoring_rate和Adjusted_contiguity在0-100之间
                    y_min = max(0, min_value * 0.5)  # 从最小值的一半开始，但不小于0
                    y_max = min(100, max_value * 1.25)  # 最大值向上留25%空间
                    if max_value < 80:
                        y_max = max_value * 1.4  # 如果最大值较小，留更多空间
            else:
                # 对于错误率指标
                y_min = 0  # 错误率从0开始
                y_max = max_value * 1.35  # 最大值向上留35%空间
                if max_value > 50:
                    y_max = max_value * 1.2  # 如果最大值很大，少留一些空间
            
            # 创建条形图 - 调整条宽度
            x_pos = np.arange(len(methods))
            bars = ax.bar(x_pos, values, color=[colors.get(m, '#999999') for m in methods], 
                         alpha=0.8, width=0.6, edgecolor='black', linewidth=1)
            
            # 设置y轴范围
            ax.set_ylim(y_min, y_max)
            
            # 添加数值标签 - 总是放在柱子顶部上方
            for bar, value in zip(bars, values):
                height = bar.get_height()
                
                # 确定标签格式
                if metric == 'Contiguity':
                    label = f'{value:.4f}'
                else:
                    label = f'{value:.2f}'
                
                # 计算标签的y位置（放在柱子顶部上方）
                # 使用相对于y轴范围的比例来确定偏移量
                y_range = y_max - y_min
                offset = y_range * 0.02  # 偏移量为y轴范围的2%
                label_y = height + offset
                
                # 确保标签不会超出图形顶部
                if label_y > y_max * 0.98:
                    label_y = height - offset  # 如果太接近顶部，放在柱子内部
                
                # 添加标签
                ax.text(bar.get_x() + bar.get_width()/2, label_y,
                       label, ha='center', va='bottom', fontsize=11, fontweight='bold',
                       color='black')
            
            # 设置图表属性
            ax.set_title(f'{metric_names[metric]}', fontsize=14, fontweight='bold', pad=15)
            
            # 设置x轴刻度和标签
            ax.set_xticks(x_pos)
            ax.set_xticklabels(methods, fontsize=12, fontweight='bold')
            
            # 添加网格线
            ax.grid(True, alpha=0.3, axis='y', linestyle='--')
            
            # 移除顶部和右侧边框
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            
            # 加粗左侧和底部边框
            ax.spines['left'].set_linewidth(1.5)
            ax.spines['bottom'].set_linewidth(1.5)
            
            # 设置y轴标签字体大小
            ax.tick_params(axis='y', labelsize=11)
            
            # 添加y轴网格线更密集
            ax.yaxis.set_major_locator(plt.MaxNLocator(8))
        
        # 设置整体标题
        main_title = f'The clustering performance of different scaffolding tools on {species_full}'
        plt.suptitle(main_title, fontsize=18, fontweight='bold', y=1.02)
        
        # 添加图例
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor=colors['ALLHiC'], alpha=0.8, edgecolor='black', linewidth=1, label='ALLHiC'),
            Patch(facecolor=colors['HapHiC'], alpha=0.8, edgecolor='black', linewidth=1, label='HapHiC'),
            Patch(facecolor=colors['ConHiC'], alpha=0.8, edgecolor='black', linewidth=1, label='ConHiC')
        ]
        fig.legend(handles=legend_elements, loc='upper center', 
                  bbox_to_anchor=(0.5, 0.98), ncol=3, fontsize=13, frameon=True,
                  framealpha=0.9, edgecolor='gray')
        
        plt.tight_layout(rect=[0, 0, 1, 0.95])
        
        # 保存图像
        filename = f'bar_chart_improved_{species_full.replace(" ", "_").replace("(", "").replace(")", "").replace("=", "")}.png'
        plt.savefig(filename, dpi=300, bbox_inches='tight', facecolor='white')
        plt.close()
        
        print(f"  已保存改进条形图: {filename}")

# 创建改进的标准化条形图
print("\n创建改进的标准化条形图:")
create_improved_bar_charts(df)

# 打印最佳值信息
print("\n" + "="*80)
print("各分组内最佳值统计:")

# 定义要评估的指标列
metric_columns = ['Contiguity', 'Inter_homo_error_rate', 'Inter_nonhomo_error_rate', 
                 'total_error_rate', 'Anchoring_rate', 'Adjusted_contiguity']

for group_idx, (group_name, group_df, _) in enumerate(groups):
    print(f"\n{group_name}:")
    group_best_values = best_values_info[group_idx][1]  # 获取最佳值信息
    
    for metric in metric_columns:
        if metric in group_best_values:
            best_idx, best_value = group_best_values[metric]
            best_row = df.loc[best_idx]
            print(f"  最佳{metric}: {best_value:.4f} ({best_row['Method']})")

# 创建详细的对比分析
print("\n" + "="*80)
print("详细对比分析:")

# 对于每个分组，详细比较
for group_idx, (group_name, group_df, _) in enumerate(groups):
    print(f"\n=== {group_name} ===")
    
    # 为每个指标找出排名（忽略NaN值）
    for metric in metric_columns:
        if metric in group_df.columns:
            # 过滤掉NaN值
            valid_data = group_df[~group_df[metric].isna()].copy()
            if len(valid_data) == 0:
                continue
            
            # 确定排序方向
            ascending = metric in ['total_error_rate', 'Inter_homo_error_rate', 'Inter_nonhomo_error_rate']
            sorted_df = valid_data.sort_values(by=metric, ascending=ascending)
            
            print(f"\n{metric}排名:")
            for i, (idx, row) in enumerate(sorted_df.iterrows()):
                rank = i + 1
                method = row['Method']
                value = row[metric]
                # 检查是否是并列第一
                is_best = (i == 0) or (sorted_df.iloc[i][metric] == sorted_df.iloc[0][metric])
                best_mark = " ★" if is_best else ""
                print(f"  {rank}. {method}: {value:.4f}{best_mark}")

# 创建综合性能排名
print("\n" + "="*80)
print("综合性能排名（每个分组内统计第一名次数）:")

# 统计每个方法在每个分组中获得第一名的次数
method_performance = {}
for group_idx, (group_name, group_df, _) in enumerate(groups):
    group_best_values = best_values_info[group_idx][1]
    
    for metric in metric_columns:
        if metric in group_best_values:
            best_idx, _ = group_best_values[metric]
            best_method = df.loc[best_idx, 'Method']
            
            if best_method not in method_performance:
                method_performance[best_method] = 0
            method_performance[best_method] += 1

# 按第一名次数排序
sorted_performance = sorted(method_performance.items(), key=lambda x: x[1], reverse=True)
print("\n方法综合表现排名（按获得第一名次数）:")
for method, count in sorted_performance:
    print(f"  {method}: {count} 次第一名")

# 保存性能排名
performance_df = pd.DataFrame(sorted_performance, columns=['Method', 'First_place_count'])
performance_csv = 'method_performance_ranking_updated.csv'
performance_df.to_csv(performance_csv, index=False, encoding='utf-8-sig')
print(f"\n方法性能排名已保存为: {performance_csv}")

# 创建按物种和方法分组的汇总表
print("\n" + "="*80)
print("各物种不同方法的平均表现:")

# 计算数值列的平均值
numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()

# 按物种和方法分组计算平均值
summary = df.groupby(['Species', 'Method'])[numeric_cols].mean().round(4)

print(summary.to_string())

# 保存汇总统计
summary_csv = 'clustering_results_summary_updated.csv'
summary.to_csv(summary_csv, encoding='utf-8-sig')
print(f"\n汇总统计已保存为: {summary_csv}")

print("\n" + "="*80)
print("所有处理完成!")
print(f"1. 原始数据已保存为: {csv_filename}")
print(f"2. 汇总表格已保存为: {summary_csv}")
print(f"3. 方法排名已保存为: {performance_csv}")
print(f"4. 主表格图像已保存为: clustering_results_table_compact.png")
print(f"5. 改进的标准化条形图已保存为: bar_chart_improved_*.png")
print("\n主要优化说明:")
print("  - 表格布局优化：增大图形宽度(22英寸)，减少列间空隙")
print("  - 列名缩写：使用更短的列名缩写(如'HomoErr%'、'NonHomoErr%')")
print("  - 字号增大：数据行字体从7号增大到10号，表头从8号增大到11号")
print("  - 数字格式：Contiguity显示4位小数，其他指标显示2位小数")
print("  - 列宽调整：自动根据内容调整列宽，减少空白区域")
print("  - 完全去除留白：通过设置表格位置[0, 0, 1, 1]和pad_inches=0去除所有边距")
print("  - 保持清晰：虽然列宽缩小，但通过增大字号保持可读性")