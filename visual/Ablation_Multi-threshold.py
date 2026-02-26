import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# 重新整理数据 - 正确的数据结构
data = {
    "Species": ["C88"] * 6 + ["SS"] * 6 + ["XJDY"] * 6,
    "Threshold": [
        "1.00", "1.03", "1.02", "1.01", "1-1.03", "1-1.02",
        "1.00", "1.03", "1.02", "1.01", "1-1.03", "1-1.02",
        "1.00", "1.01", "1.02", "1.03", "1-1.03", "1-1.02"
    ],
    "Syntenic_contigs": [57.2619, 55.9716, 56.6696, 57.5698, 57.5698, 57.5698,
                         60.3874, 61.8814, 63.8515, 63.4759, 63.8523, 63.344,
                         61.3037, 61.229, 60.735, 59.2787, 61.0336, 61.1671],
    "Unanchored_contigs": [8.8432, 8.8432, 8.8432, 8.8432, 8.8432, 8.8432,
                           9.342, 9.342, 9.342, 9.342, 9.342, 9.342,
                           11.2683, 11.2683, 11.2683, 11.2683, 11.2683, 11.2683],
    "Newly_anchored_contigs": [1.1963, 1.1963, 1.1963, 1.1963, 1.1963, 1.1963,
                               0.0271, 0.0271, 0.0271, 0.0271, 0.0271, 0.0271,
                               8.0157, 8.0157, 8.0157, 8.0157, 8.0157, 8.0157],
    "Translocation_contigs": [9.404, 9.4049, 9.4049, 9.4049, 9.4049, 9.4049,
                              8.7501, 8.7501, 8.7501, 8.7501, 8.7501, 8.7501,
                              10.5326, 10.5326, 10.5326, 10.5326, 10.5326, 10.5326],
    "Relocation_contigs": [9.929, 11.3179, 10.2933, 9.9135, 9.9831, 9.9135,
                           9.9795, 7.3462, 4.8305, 5.9495, 6.0907, 6.2191,
                           1.7625, 1.635, 1.9029, 2.707, 1.7362, 1.7386],
    "Inversion_contigs": [8.71648, 8.6025, 8.3072, 8.3121, 8.2425, 8.3121,
                          7.464, 8.5846, 8.5609, 8.1664, 8.0251, 8.1341,
                          4.3783, 4.4356, 4.5193, 4.717, 4.4873, 4.4507],
    "Inversion_relocation_contigs": [4.6479, 4.6632, 5.2849, 4.7598, 4.7598, 4.7598,
                                     4.0494, 4.0683, 4.6375, 4.2886, 3.9122, 4.1832,
                                     2.7386, 2.8835, 3.0259, 3.4804, 2.926, 2.8268]
}

# 创建DataFrame
df = pd.DataFrame(data)

# 修改表格展示（三个物种分开）
fig, axes = plt.subplots(3, 1, figsize=(22, 18))
fig.suptitle('Genome Assembly Quality - Complete Results Table for Three Species', 
              fontsize=24, fontweight='bold', y=0.98)

# 定义表格列名
column_names = ['Threshold', 'Syntenic\n(%)', 'Unanchored\n(%)', 
                'Newly\nAnchored\n(%)', 'Translocation\n(%)', 
                'Relocation\n(%)', 'Inversion\n(%)', 'Inv+Reloc\n(%)']

# 为每个物种创建一个表格
for idx, species in enumerate(['C88', 'SS', 'XJDY']):
    ax = axes[idx]
    ax.axis('tight')
    ax.axis('off')
    
    # 获取该物种的数据
    species_data = df[df['Species'] == species].copy()
    
    # 按阈值顺序排序
    threshold_order_specific = ['1.00', '1.01', '1.02', '1.03', '1-1.02', '1-1.03']
    ordered_data = []
    
    for threshold in threshold_order_specific:
        threshold_data = species_data[species_data['Threshold'] == threshold]
        if not threshold_data.empty:
            row = threshold_data.iloc[0]
            ordered_data.append([
                threshold,
                row['Syntenic_contigs'],  # 保存原始数值用于比较
                row['Unanchored_contigs'],
                row['Newly_anchored_contigs'],
                row['Translocation_contigs'],
                row['Relocation_contigs'],
                row['Inversion_contigs'],
                row['Inversion_relocation_contigs']
            ])
        else:
            # 对于缺失的阈值，使用最接近的值
            if threshold == '1.01' and species == 'C88':
                # C88没有1.01，使用1.00的数据
                row = species_data[species_data['Threshold'] == '1.00'].iloc[0]
                ordered_data.append([
                    threshold + '*',
                    row['Syntenic_contigs'],
                    row['Unanchored_contigs'],
                    row['Newly_anchored_contigs'],
                    row['Translocation_contigs'],
                    row['Relocation_contigs'],
                    row['Inversion_contigs'],
                    row['Inversion_relocation_contigs']
                ])
    
    # 找到每列的最优值（对于不同的指标，最优的定义不同）
    # Syntenic: 越高越好
    # Unanchored: 越低越好
    # Newly_anchored: 对于SS越低越好，对于其他可能不同
    # Translocation, Relocation, Inversion, Inv+Reloc: 越低越好
    
    # 获取每列的最优值索引
    best_indices = {}
    
    # Syntenic列：最高值最优
    syntenic_values = [row[1] for row in ordered_data]
    best_syntenic_idx = np.argmax(syntenic_values)
    best_indices[1] = best_syntenic_idx
    
    # Unanchored列：最低值最优
    unanchored_values = [row[2] for row in ordered_data]
    best_unanchored_idx = np.argmin(unanchored_values)
    best_indices[2] = best_unanchored_idx
    
    # Newly_anchored列：对于SS，最低值最优；对于其他，根据情况判断
    if species == 'SS':
        newly_values = [row[3] for row in ordered_data]
        best_newly_idx = np.argmin(newly_values)
    else:
        # 对于C88和XJDY，暂时也按最低值处理
        newly_values = [row[3] for row in ordered_data]
        best_newly_idx = np.argmin(newly_values)
    best_indices[3] = best_newly_idx
    
    # 错误类型列：最低值最优
    error_columns = [4, 5, 6, 7]  # Translocation, Relocation, Inversion, Inv+Reloc
    for col in error_columns:
        values = [row[col] for row in ordered_data]
        best_idx = np.argmin(values)
        best_indices[col] = best_idx
    
    # 将数值转换为字符串格式，并标记最优值
    formatted_data = []
    for i, row in enumerate(ordered_data):
        formatted_row = []
        for j, value in enumerate(row):
            if j == 0:  # 阈值列
                formatted_row.append(str(value))
            else:  # 数值列
                # 检查是否为最优值
                is_best = (i == best_indices.get(j, -1))
                if is_best:
                    formatted_row.append(f"{value:.4f}")  # 最优值将用加粗显示
                else:
                    formatted_row.append(f"{value:.4f}")
        formatted_data.append(formatted_row)
    
    # 创建表格
    table = ax.table(cellText=formatted_data,
                     colLabels=column_names,
                     cellLoc='center',
                     loc='center',
                     colWidths=[0.1] + [0.1125]*7)
    
    table.auto_set_font_size(False)
    table.set_fontsize(12)
    table.scale(1.3, 2.2)
    
    # 设置标题行颜色（深蓝色）
    for j in range(8):
        table[(0, j)].set_facecolor('#2c3e50')
        table[(0, j)].set_text_props(color='white', weight='bold', fontsize=13)
    
    # 设置数据行样式
    for i in range(1, len(formatted_data)+1):
        for j in range(8):
            # 所有单元格背景色为白色
            table[(i, j)].set_facecolor('white')
            
            # 所有文本颜色为黑色
            table[(i, j)].set_text_props(color='black', fontsize=12)
            
            # 如果是数值列且为最优值，则加粗显示
            if j >= 1 and j <= 7:
                # 检查是否为最优值
                original_i = i - 1  # 转换为原始数据索引
                if original_i == best_indices.get(j, -1):
                    table[(i, j)].set_text_props(color='black', weight='bold', fontsize=12)
    
    # 设置物种标题（黑色）
    species_title = f'{species} Assembly Quality Results'
    if species == 'C88':
        species_title += ' (N50=2M)'
    
    ax.set_title(species_title, fontsize=16, fontweight='bold', pad=30, 
                 color='black')

plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig('Assembly_Quality_Complete_Table_Modified.png', dpi=300, bbox_inches='tight')
plt.show()

print("修改后的表格图已保存为: Assembly_Quality_Complete_Table_Modified.png")

# 创建一个更简洁的版本，使用交替行颜色但最后一行为白色
fig2, axes2 = plt.subplots(3, 1, figsize=(22, 18))
fig2.suptitle('Genome Assembly Quality - Complete Results Table for Three Species', 
              fontsize=24, fontweight='bold', y=0.98)

# 为每个物种创建一个表格（交替行颜色版本）
for idx, species in enumerate(['C88', 'SS', 'XJDY']):
    ax = axes2[idx]
    ax.axis('tight')
    ax.axis('off')
    
    # 获取该物种的数据
    species_data = df[df['Species'] == species].copy()
    
    # 按阈值顺序排序
    threshold_order_specific = ['1.00', '1.01', '1.02', '1.03', '1-1.02', '1-1.03']
    ordered_data = []
    
    for threshold in threshold_order_specific:
        threshold_data = species_data[species_data['Threshold'] == threshold]
        if not threshold_data.empty:
            row = threshold_data.iloc[0]
            ordered_data.append([
                threshold,
                row['Syntenic_contigs'],  # 保存原始数值用于比较
                row['Unanchored_contigs'],
                row['Newly_anchored_contigs'],
                row['Translocation_contigs'],
                row['Relocation_contigs'],
                row['Inversion_contigs'],
                row['Inversion_relocation_contigs']
            ])
    
    # 找到每列的最优值
    best_indices = {}
    
    # Syntenic列：最高值最优
    syntenic_values = [row[1] for row in ordered_data]
    best_syntenic_idx = np.argmax(syntenic_values)
    best_indices[1] = best_syntenic_idx
    
    # Unanchored列：最低值最优
    unanchored_values = [row[2] for row in ordered_data]
    best_unanchored_idx = np.argmin(unanchored_values)
    best_indices[2] = best_unanchored_idx
    
    # Newly_anchored列：最低值最优
    newly_values = [row[3] for row in ordered_data]
    best_newly_idx = np.argmin(newly_values)
    best_indices[3] = best_newly_idx
    
    # 错误类型列：最低值最优
    error_columns = [4, 5, 6, 7]  # Translocation, Relocation, Inversion, Inv+Reloc
    for col in error_columns:
        values = [row[col] for row in ordered_data]
        best_idx = np.argmin(values)
        best_indices[col] = best_idx
    
    # 将数值转换为字符串格式
    formatted_data = []
    for i, row in enumerate(ordered_data):
        formatted_row = []
        for j, value in enumerate(row):
            if j == 0:  # 阈值列
                formatted_row.append(str(value))
            else:  # 数值列
                formatted_row.append(f"{value:.4f}")
        formatted_data.append(formatted_row)
    
    # 创建表格
    table = ax.table(cellText=formatted_data,
                     colLabels=column_names,
                     cellLoc='center',
                     loc='center',
                     colWidths=[0.1] + [0.1125]*7)
    
    table.auto_set_font_size(False)
    table.set_fontsize(12)
    table.scale(1.3, 2.2)
    
    # 设置标题行颜色（深蓝色）
    for j in range(8):
        table[(0, j)].set_facecolor('#2c3e50')
        table[(0, j)].set_text_props(color='white', weight='bold', fontsize=13)
    
    # 设置数据行样式 - 使用交替行颜色，但最后一行保持白色
    for i in range(1, len(formatted_data)+1):
        for j in range(8):
            # 如果是最后一行，背景为白色
            if i == len(formatted_data):
                table[(i, j)].set_facecolor('white')
            else:
                # 其他行使用交替颜色
                if i % 2 == 0:
                    table[(i, j)].set_facecolor('#f8f9fa')  # 浅灰色
                else:
                    table[(i, j)].set_facecolor('white')   # 白色
            
            # 所有文本颜色为黑色
            table[(i, j)].set_text_props(color='black', fontsize=12)
            
            # 如果是数值列且为最优值，则加粗显示
            if j >= 1 and j <= 7:
                # 检查是否为最优值
                original_i = i - 1  # 转换为原始数据索引
                if original_i == best_indices.get(j, -1):
                    table[(i, j)].set_text_props(color='black', weight='bold', fontsize=12)
    
    # 设置物种标题（黑色）
    species_title = f'{species} Assembly Quality Results'
    if species == 'C88':
        species_title += ' (N50=2M)'
    
    ax.set_title(species_title, fontsize=16, fontweight='bold', pad=30, 
                 color='black')

plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig('Assembly_Quality_Complete_Table_Alternate.png', dpi=300, bbox_inches='tight')
plt.show()

print("交替行颜色版本的表格图已保存为: Assembly_Quality_Complete_Table_Alternate.png")

print("\n" + "="*70)
print("表格修改说明:")
print("="*70)
print("1. 所有数值使用黑色字体")
print("2. 每列的最优值用加粗显示")
print("3. 表格标题使用黑色")
print("4. 生成了两个版本:")
print("   - 版本1: 所有数据行为白色背景")
print("   - 版本2: 交替行颜色，但最后一行保持白色")
print("="*70)