import numpy as np
import matplotlib.pyplot as plt
import os

# 数据准备 - 包含所有4个工具和所有数据集
# 使用 None 表示缺失数据
data = {
    'ALLHiC': {
        'C88(2M)': [46.0673, 16.3802, 1.4047, 19.2147, 7.6956, 6.4205, 2.8166],
        'C88(0.8M)': [43.0654, 20.8614, 1.2506, 16.72, 8.7514, 5.0462, 4.3046],
        'SS': [59.7134, 0.1663, 1.1653, 23.4162, 3.4706, 7.2821, 4.7859],
        'XJDY': [19.5894, 7.7619, 9.9995, 56.598, 1.6375, 1.5206, 2.8928],
        'ZM-4': [26.4088, 1.1927, 6.3062, 38.6665, 9.2421, 8.194, 9.9892],
        'cs': [11.7541, 19.1868, 0.2432, 42.65, 9.0222, 2.8755, 14.2679],
        'ap': None
    },
    'HapHiC': {
        'C88(2M)': [55.5551, 10.252, 1.2111, 12.4265, 8.0688, 7.7714, 4.7146],
        'C88(0.8M)': [51.8432, 15.5156, 0.9271, 11.6222, 8.7582, 5.2423, 6.0908],
        'SS': [60.9713, 9.3076, 0.0936, 12.5326, 5.3461, 6.9315, 4.817],
        'XJDY': [61.2931, 11.288, 7.9993, 10.5226, 1.7656, 4.3877, 2.7433],
        'ZM-4': [42.1392, 5.1487, 5.3584, 16.1181, 9.4798, 11.3425, 10.4129],
        'cs': [31.1367, 21.4091, 0.001, 35.6152, 3.1962, 4.0061, 4.6355],
        'ap': [34.8641, 6.2913, 4.1401, 15.7736, 11.9817, 14.4828, 12.466]
    },
    'ConHiC (cluster)': {
        'C88(2M)': [57.2619, 8.8432, 1.1963, 9.404, 9.929, 8.71648, 4.6479],
        'C88(0.8M)': [54.0403, 13.0403, 1.0127, 10.9585, 9.0885, 5.5918, 6.2675],
        'SS': [60.3874, 9.342, 0.0271, 8.7501, 9.9795, 7.464, 4.0494],
        'XJDY': [61.3037, 11.2683, 8.0157, 10.5326, 1.7625, 4.3783, 2.7386],
        'ZM-4': [42.0324, 4.9225, 5.3997, 15.3528, 9.8666, 11.6961, 10.7296],
        'cs': [31.2574, 21.3717, 0.001, 35.608, 3.2087, 3.9447, 4.6082],
        'ap': [35.6971, 5.6955, 4.2409, 15.9904, 11.3626, 15.0478, 11.9653]
    },
    'ConHiC (cluster+sort)': {
        'C88(2M)': [57.5698, 8.8432, 1.1963, 9.4049, 9.9831, 8.2425, 4.7598],
        'C88(0.8M)': [53.3584, 13.0403, 1.0127, 10.9585, 9.4122, 6.302, 5.9156],
        'SS': [63.8523, 9.342, 0.0271, 8.7501, 6.0907, 8.0251, 3.9122],
        'XJDY': [61.0336, 11.2683, 8.0157, 10.5326, 1.7362, 4.4873, 2.926],
        'ZM-4': [42.5891, 4.9225, 5.3997, 15.3528, 9.5636, 11.6863, 10.4857],
        'cs': [31.2832, 21.3717, 0.001, 35.608, 3.224, 3.8984, 4.6134],
        'ap': [35.8503, 5.6955, 4.2409, 15.9904, 11.3654, 14.9706, 11.8866]
    }
}

# 指标名称和优化方向
categories = [
    'Syntenic contigs\n(higher is better)',
    'Unanchored contigs\n(lower is better)', 
    'Newly anchored contigs\n(higher is better)',
    'Translocation contigs\n(lower is better)',
    'Relocation contigs\n(lower is better)',
    'Inversion contigs\n(lower is better)',
    'Inversion and relocation\ncontigs (lower is better)'
]

# 简化的指标名称用于表格
short_categories = [
    'Syntenic\n(higher ↑)',
    'Unanchored\n(lower ↓)', 
    'Newly anchored\n(higher ↑)',
    'Translocation\n(lower ↓)',
    'Relocation\n(lower ↓)',
    'Inversion\n(lower ↓)',
    'Inversion+Reloc\n(lower ↓)'
]

# 优化方向：1表示越大越好，-1表示越小越好
optimization_direction = [1, -1, 1, -1, -1, -1, -1]

N = len(categories)

# 工具配置 - 4个工具
tools = ['ALLHiC', 'HapHiC', 'ConHiC (cluster)', 'ConHiC (cluster+sort)']

# 所有数据集
datasets = ['C88(2M)', 'C88(0.8M)', 'SS', 'XJDY', 'ZM-4', 'cs', 'ap']

# 创建保存目录
output_dir = 'performance_tables_all_datasets'
os.makedirs(output_dir, exist_ok=True)

# 设置全局样式
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 9

print("正在分析数据并确定最佳值...")

# 首先确定每个数据集每个指标的最佳值
# 只比较有值的工具，如果有多个最佳值相等，只选择最靠上的（按tools列表顺序）
best_values = {}  # 格式: {dataset: {metric_index: (best_value, tool_index)}}

for dataset_idx, dataset in enumerate(datasets):
    best_values[dataset] = {}
    for metric_idx in range(N):
        # 收集所有工具在该指标上的有效值
        values = []
        tool_indices = []
        for tool_idx, tool in enumerate(tools):
            tool_data = data[tool][dataset]
            if tool_data is not None:
                value = tool_data[metric_idx]
                values.append(value)
                tool_indices.append(tool_idx)
        
        # 如果有有效值，确定最佳值
        if values:
            if optimization_direction[metric_idx] == 1:  # 越大越好
                best_value = max(values)
                # 找到第一个达到最佳值的工具（最靠上）
                for i, (val, tool_idx) in enumerate(zip(values, tool_indices)):
                    if abs(val - best_value) < 0.0001:  # 考虑浮点误差
                        best_tool_idx = tool_idx
                        break
            else:  # 越小越好
                best_value = min(values)
                # 找到第一个达到最佳值的工具（最靠上）
                for i, (val, tool_idx) in enumerate(zip(values, tool_indices)):
                    if abs(val - best_value) < 0.0001:  # 考虑浮点误差
                        best_tool_idx = tool_idx
                        break
            
            best_values[dataset][metric_idx] = (best_value, best_tool_idx)

print("\n正在生成数据表格（最佳值加粗，不同数据集间有间隔）...")
fig_table, ax_table = plt.subplots(figsize=(16, 14))
ax_table.axis('tight')
ax_table.axis('off')

# 准备表格数据，在不同数据集之间添加空行
# 去掉ap;allhic这一行
table_data = []
headers = ['Dataset', 'Tool'] + short_categories

# 存储哪些单元格需要加粗
bold_cells = []  # 格式: [(row_idx, col_idx), ...]

# 为每个数据集和工具添加数据，并标记最佳值
# 同时在每个数据集之后添加一个空行作为间隔
row_counter = 0
for dataset_idx, dataset in enumerate(datasets):
    # 为当前数据集添加所有工具的数据
    dataset_has_any_data = False
    
    for tool_idx, tool in enumerate(tools):
        # 去掉ap;allhic这一行
        if dataset == 'ap' and tool == 'ALLHiC':
            continue
            
        tool_data = data[tool][dataset]
        formatted_values = []
        
        # 检查该工具在该数据集上是否有任何有效数据
        if tool_data is None:
            dataset_has_any_data = True
            # 全部用"-"填充
            formatted_values = ['-'] * N
            row = [dataset, tool] + formatted_values
            table_data.append(row)
            row_counter += 1
        else:
            dataset_has_any_data = True
            # 格式化每个值，并标记最佳值
            for metric_idx, value in enumerate(tool_data):
                value_str = f'{value:.3f}'
                formatted_values.append(value_str)
                
                # 检查是否是当前数据集当前指标的最佳值
                if (dataset in best_values and metric_idx in best_values[dataset]):
                    best_value, best_tool_idx = best_values[dataset][metric_idx]
                    if tool_idx == best_tool_idx and abs(value - best_value) < 0.0001:
                        # 记录这个单元格需要加粗
                        bold_cells.append((row_counter + 1, metric_idx + 2))  # +1因为表头行，+2因为Dataset和Tool列
            
            row = [dataset, tool] + formatted_values
            table_data.append(row)
            row_counter += 1
    
    # 如果数据集有任何数据，在当前数据集的所有工具数据后添加一个空行（除了最后一个数据集）
    if dataset_has_any_data and dataset_idx < len(datasets) - 1:
        # 创建空行
        empty_row = [''] * len(headers)
        table_data.append(empty_row)
        row_counter += 1

print(f"\n需要加粗的单元格数量: {len(bold_cells)}")
print("需要加粗的单元格位置（行,列）示例:", bold_cells[:10])  # 只显示前10个

# 创建表格
table = ax_table.table(cellText=table_data, colLabels=headers, loc='center', cellLoc='center')
table.auto_set_font_size(False)
table.set_fontsize(7.5)
table.scale(1, 1.8)

# 设置标题行样式
for j in range(len(headers)):
    table[(0, j)].set_facecolor('#404040')
    table[(0, j)].get_text().set_color('white')
    table[(0, j)].get_text().set_weight('bold')

# 设置交替行颜色和工具颜色
tool_color_map = {
    'ALLHiC': '#e6f2ff',        # 浅蓝色
    'HapHiC': '#e6ffe6',        # 浅绿色
    'ConHiC (cluster)': '#fff0e6',  # 浅橙色
    'ConHiC (cluster+sort)': '#ffe6e6'  # 浅红色
}

# 为每个单元格设置样式
current_row = 1  # 从第1行开始（0是标题行）
table_row_idx = 0  # table_data的行索引

for i, row_data in enumerate(table_data):
    if i >= len(table_data):
        break
    
    # 跳过空行
    if not row_data[0]:  # 如果是空行（第一列为空）
        # 设置空行的样式（浅灰色背景）
        for j in range(len(headers)):
            table[(current_row, j)].set_facecolor('#f8f8f8')
            table[(current_row, j)].get_text().set_color('#f8f8f8')
        current_row += 1
        table_row_idx += 1
        continue
    
    dataset = row_data[0]
    tool = row_data[1]
    
    color = tool_color_map.get(tool, '#ffffff')
    
    # 为整个行设置背景色
    for j in range(len(headers)):
        table[(current_row, j)].set_facecolor(color)
    
    # 为数据集列添加灰色背景区分（每个数据集的第一个工具）
    is_first_tool_in_dataset = (i == 0 or not table_data[i-1][0] or table_data[i-1][0] != dataset)
    if is_first_tool_in_dataset:
        table[(current_row, 0)].set_facecolor('#f0f0f0')
        table[(current_row, 0)].get_text().set_weight('bold')
    
    # 设置数据单元格样式（先设置默认样式）
    for metric_idx in range(N):
        col_idx = metric_idx + 2  # 前两列是Dataset和Tool
        cell_value = row_data[col_idx]
        
        if cell_value == '-':
            # 对于"-"值，设置为灰色斜体
            table[(current_row, col_idx)].get_text().set_color('#666666')
            table[(current_row, col_idx)].get_text().set_style('italic')
        else:
            # 默认样式
            table[(current_row, col_idx)].get_text().set_color('#000000')
    
    current_row += 1
    table_row_idx += 1

# 设置最佳值的加粗样式
print("\n正在设置最佳值的加粗样式...")
bold_count = 0
for (row_idx, col_idx) in bold_cells:
    try:
        # 直接设置单元格文本的加粗样式
        table[(row_idx, col_idx)].get_text().set_weight('bold')
        table[(row_idx, col_idx)].get_text().set_color('#000000')
        bold_count += 1
    except Exception as e:
        print(f"  警告: 无法加粗行{row_idx}列{col_idx}: {e}")

print(f"成功加粗了 {bold_count} 个单元格")

# 添加边框
for i in range(len(table_data) + 1):
    for j in range(len(headers)):
        table[(i, j)].set_edgecolor('#cccccc')

# 优化列宽
table.auto_set_column_width([i for i in range(len(headers))])

#ax_table.set_title('Detailed Performance Data (Original Values)\nBold values indicate best performance for each metric\n"-" indicates data not available (ALLHiC cannot measure on ap; other tools not tested yet)\n↑ higher is better, ↓ lower is better', 
#                  fontsize=12, weight='bold', y=1.02)
plt.tight_layout()

# 保存表格
output_path_table = os.path.join(output_dir, 'performance_data_table_4tools_bold_with_spacing.png')
plt.savefig(output_path_table, dpi=300, bbox_inches='tight', facecolor='white')
print(f"\n数据表格已保存至: {output_path_table}")

# 同时保存为CSV格式以便进一步处理
print("\n正在生成CSV格式的数据表格...")
csv_file_path = os.path.join(output_dir, 'performance_data_table_4tools.csv')

# 准备CSV数据
csv_data = []
# 添加表头
csv_header = ['Dataset', 'Tool'] + short_categories
csv_data.append(csv_header)

# 添加数据行
for i, row_data in enumerate(table_data):
    if row_data[0]:  # 如果不是空行
        csv_data.append(row_data)

# 写入CSV文件
with open(csv_file_path, 'w', encoding='utf-8') as f:
    for row in csv_data:
        f.write(','.join([str(x) for x in row]) + '\n')

print(f"CSV表格已保存至: {csv_file_path}")

plt.show()

# 生成详细的最佳值统计
print(f"\n{'='*70}")
print("详细的最佳值分布统计:")
print(f"{'='*70}")

# 按数据集统计
for dataset in ['C88(2M)', 'C88(0.8M)', 'SS', 'XJDY', 'ZM-4', 'cs', 'ap']:
    if dataset in best_values and best_values[dataset]:
        print(f"\n{dataset}:")
        for metric_idx in range(N):
            if metric_idx in best_values[dataset]:
                best_value, best_tool_idx = best_values[dataset][metric_idx]
                tool_name = tools[best_tool_idx]
                metric_name = short_categories[metric_idx]
                print(f"  {metric_name:<20} {best_value:.4f} ({tool_name})")

# 按工具统计
print(f"\n{'='*70}")
print("最佳值总数统计:")
print(f"{'='*70}")
tool_best_counts = {tool: 0 for tool in tools}
total_comparisons = 0

for dataset in datasets:
    for metric_idx in range(N):
        if dataset in best_values and metric_idx in best_values[dataset]:
            best_value, best_tool_idx = best_values[dataset][metric_idx]
            tool_best_counts[tools[best_tool_idx]] += 1
            total_comparisons += 1

print(f"总比较次数: {total_comparisons}")
print(f"{'-'*40}")
for tool in tools:
    count = tool_best_counts[tool]
    percentage = (count / total_comparisons * 100) if total_comparisons > 0 else 0
    print(f"{tool:<25} {count:>3} ({percentage:>5.1f}%)")

print(f"\n{'='*70}")
print("所有表格已成功保存！")
print(f"保存目录: {output_dir}")
print(f"\n输出文件:")
print(f"  - {output_path_table} (最佳值加粗，不同数据集间有间隔)")
print(f"  - {csv_file_path} (CSV格式)")
print(f"{'='*70}")