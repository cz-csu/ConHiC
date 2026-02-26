import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

# 数据
data = {
    "Resolution": ["0.04M", "0.05M", "0.1M", "0.25M", "0.5M", "0.75M", "1M", "1.05M", "1.06M", "1.25M", "1.5M", "2M", "2.5M", "3M"],
    "Contiguity": [0.8428, 0.8740, 0.7986, 0.8528, 0.8588, 0.8652, 0.8697, 0.8711, 0.8635, 0.8678, 0.8635, 0.8472, 0.8033, 0.8125],
    "Inter_homo_error_rate": [5.1087, 7.6874, 10.2965, 5.7521, 6.6713, 6.0293, 5.5189, 5.6960, 6.1586, 5.7235, 5.9363, 6.0354, 6.1617, 6.6769],
    "Inter_nonhomo_error_rate": [12.6540, 4.7045, 11.8316, 6.1260, 5.0437, 5.1973, 5.1750, 5.0247, 5.1193, 5.0819, 6.8794, 5.4297, 5.7005, 7.8670],
    "total_error_rate": [17.7627, 12.3919, 22.1281, 11.8781, 11.7150, 11.2266, 10.6939, 10.7207, 11.2779, 10.8054, 12.8157, 11.4651, 11.8622, 14.5439],
    "Anchoring rate": [88.2720, 89.5491, 77.8720, 89.3761, 89.8024, 89.7129, 89.8151, 89.8722, 89.8770, 89.4552, 89.6511, 88.5367, 84.4834, 85.8544],
    "Adjusted contiguity": [74.3956, 78.2659, 62.1886, 76.2199, 77.1223, 77.6196, 78.1122, 78.2877, 77.6088, 77.6292, 77.4137, 75.0083, 67.8655, 69.7567]
}

multi_data = {
    "Method": [
        "1M & 1.05M & 1.06M CI",
        "1M & 1.25M CI",
        "0.5M & 1M & 2M CI",
        "0.5M & 1M & 2M MC",
        "0.75M & 1M & 1.25M CI",
        "0.75M & 1M & 1.25M MC",
        "0.5M & 1M & 1.5M & 2M CI",
        "0.5M & 1M & 1.5M & 2M MC"
    ],
    "Contiguity": [0.8663, 0.8680, 0.8719, 0.8755, 0.8693, 0.8673, 0.8731, 0.8641],
    "Inter_homo_error_rate": [5.8795, 5.1815, 4.8830, 4.9074, 5.0675, 5.8257, 4.8267, 5.6510],
    "Inter_nonhomo_error_rate": [5.1454, 4.8269, 6.4721, 6.7378, 4.7609, 5.1467, 6.4108, 7.0209],
    "total_error_rate": [11.0249, 10.0084, 11.3551, 11.6452, 9.8284, 10.9724, 11.2375, 12.6719],
    "Anchoring rate": [89.9428, 88.4299, 88.7840, 89.0381, 87.8290, 89.9192, 88.6663, 90.0025],
    "Adjusted contiguity": [77.9174, 76.7572, 77.4108, 77.9529, 76.3497, 77.9869, 77.4145, 77.7712]
}

# 创建DataFrame
df_single = pd.DataFrame(data)
df_multi = pd.DataFrame(multi_data)

# 1. 创建折线图
fig1, axes = plt.subplots(2, 2, figsize=(18, 12))
fig1.suptitle('C88 N50=2M - Performance Metrics Across Resolutions', fontsize=16, fontweight='bold', y=0.98)

# 准备x轴标签（单分辨率 + 多分辨率）
single_labels = df_single['Resolution'].tolist()
multi_labels = df_multi['Method'].tolist()
all_labels = single_labels + multi_labels
x_positions = np.arange(len(all_labels))

# 定义指标和对应的颜色
metrics = [
    ('Adjusted contiguity', 'Adj Contiguity', 'blue', 'o-'),
    ('Anchoring rate', 'Anchoring Rate', 'green', 's-'),
    ('Contiguity', 'Contiguity', 'red', '^-'),
    ('total_error_rate', 'Total Error Rate', 'orange', 'd-')
]

# 准备数据（标准化到相似范围以便比较）
for idx, (ax, (metric_key, metric_name, color, marker)) in enumerate(zip(axes.flatten(), metrics)):
    # 单分辨率数据
    single_values = df_single[metric_key].values
    if metric_key == 'Contiguity':  # Contiguity需要乘以100
        single_values = single_values * 100
    
    # 多分辨率数据
    multi_values = df_multi[metric_key].values
    if metric_key == 'Contiguity':  # Contiguity需要乘以100
        multi_values = multi_values * 100
    
    # 合并数据
    all_values = np.concatenate([single_values, multi_values])
    
    # 绘制折线
    ax.plot(x_positions, all_values, marker, color=color, linewidth=2, markersize=6, label=metric_name)
    
    # 添加垂直线区分单分辨率和多分辨率
    ax.axvline(x=len(single_labels)-0.5, color='gray', linestyle='--', alpha=0.7, linewidth=1)
    
    # 设置x轴标签
    ax.set_xticks(x_positions)
    ax.set_xticklabels(all_labels, rotation=45, ha='right', fontsize=9)
    
    # 设置y轴标签
    if metric_key == 'Contiguity':
        ax.set_ylabel('Contiguity (%)', fontsize=10)
    elif metric_key == 'total_error_rate':
        ax.set_ylabel('Error Rate (%)', fontsize=10)
    else:
        ax.set_ylabel(f'{metric_name} (%)', fontsize=10)
    
    # 设置标题
    ax.set_title(f'{metric_name} Comparison', fontsize=12, fontweight='bold')
    
    # 添加网格
    ax.grid(True, alpha=0.3)
    
    # 添加图例
    ax.legend(loc='best')
    
    # 设置y轴范围
    if metric_key == 'Adjusted contiguity':
        ax.set_ylim(60, 85)
    elif metric_key == 'Anchoring rate':
        ax.set_ylim(75, 95)
    elif metric_key == 'Contiguity':
        ax.set_ylim(75, 95)
    elif metric_key == 'total_error_rate':
        ax.set_ylim(5, 25)

# 调整布局
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig('C88_Performance_Line_Charts.png', dpi=300, bbox_inches='tight')
plt.show()

print(f"折线图已保存为: C88_Performance_Line_Charts.png")

# 2. 修正C88_All_Results_Table.png表格的标红问题
fig2, ax = plt.subplots(figsize=(20, 12))
ax.axis('tight')
ax.axis('off')

# 准备完整数据
all_data = []

# 添加单分辨率标题
all_data.append(["--- SINGLE RESOLUTION ---", "", "", "", "", "", ""])

# 添加单分辨率列标题（这行应该标红）
all_data.append(["Resolution", "Contiguity", "Homo Error", "Non-homo Error", 
                 "Total Error", "Anchoring Rate", "Adj Contiguity"])

# 添加单分辨率数据
for _, row in df_single.iterrows():
    all_data.append([
        row['Resolution'],
        f"{row['Contiguity']:.4f}",
        f"{row['Inter_homo_error_rate']:.4f}",
        f"{row['Inter_nonhomo_error_rate']:.4f}",
        f"{row['total_error_rate']:.4f}",
        f"{row['Anchoring rate']:.4f}",
        f"{row['Adjusted contiguity']:.4f}"
    ])

# 添加分隔行
all_data.append(["", "", "", "", "", "", ""])
all_data.append(["--- MULTI RESOLUTION ---", "", "", "", "", "", ""])

# 添加多分辨率列标题（这行应该标红）
all_data.append(["Method", "Contiguity", "Homo Error", "Non-homo Error", 
                 "Total Error", "Anchoring Rate", "Adj Contiguity"])

# 添加多分辨率数据
for _, row in df_multi.iterrows():
    all_data.append([
        row['Method'],
        f"{row['Contiguity']:.4f}",
        f"{row['Inter_homo_error_rate']:.4f}",
        f"{row['Inter_nonhomo_error_rate']:.4f}",
        f"{row['total_error_rate']:.4f}",
        f"{row['Anchoring rate']:.4f}",
        f"{row['Adjusted contiguity']:.4f}"
    ])

# 创建完整表格
full_table = ax.table(cellText=all_data,
                      cellLoc='center',
                      loc='center')

# 设置表格样式
full_table.auto_set_font_size(False)
full_table.set_fontsize(9)
full_table.scale(1.3, 1.5)

# 正确设置颜色：列标题行（第1行和第len(df_single)+5行）标红
for i in range(len(all_data)):
    if "SINGLE RESOLUTION" in str(all_data[i][0]) or "MULTI RESOLUTION" in str(all_data[i][0]):
        # 分区标题行 - 灰色
        for j in range(7):
            full_table[(i, j)].set_facecolor('#6c757d')
            full_table[(i, j)].set_text_props(color='white', weight='bold', fontsize=11)
    elif i == 1:  # 单分辨率列标题行 - 蓝色
        for j in range(7):
            full_table[(i, j)].set_facecolor('#4a86e8')
            full_table[(i, j)].set_text_props(color='white', weight='bold', fontsize=10)
    elif i == len(df_single) + 5:  # 多分辨率列标题行 - 红色
        for j in range(7):
            full_table[(i, j)].set_facecolor('#d9534f')
            full_table[(i, j)].set_text_props(color='white', weight='bold', fontsize=10)
    else:  # 数据行
        for j in range(7):
            if i > 1 and i < len(df_single) + 2:  # 单分辨率数据行
                if i % 2 == 0:
                    full_table[(i, j)].set_facecolor('#f8f9fa')
            elif i > len(df_single) + 5:  # 多分辨率数据行
                if (i - len(df_single) - 6) % 2 == 0:
                    full_table[(i, j)].set_facecolor('#f8f9fa')

plt.title('C88 N50=2M - All Results in One Table', fontsize=16, fontweight='bold', pad=20)
plt.tight_layout()

full_filename = "C88_All_Results_Table_Fixed.png"
plt.savefig(full_filename, dpi=300, bbox_inches='tight')
plt.show()

print(f"修正后的完整表格已保存为: {full_filename}")

# 3. 修正C88_Compact_Results_Table.png，将多尺度表格往上挪
fig3, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 12))

# 单分辨率表格
ax1.axis('tight')
ax1.axis('off')

# 准备单分辨率数据
single_table_data = []
for _, row in df_single.iterrows():
    single_table_data.append([
        row['Resolution'],
        f"{row['Contiguity']:.4f}",
        f"{row['Inter_homo_error_rate']:.4f}",
        f"{row['Inter_nonhomo_error_rate']:.4f}",
        f"{row['total_error_rate']:.4f}",
        f"{row['Anchoring rate']:.4f}",
        f"{row['Adjusted contiguity']:.4f}"
    ])

single_table = ax1.table(cellText=single_table_data,
                        colLabels=['Res', 'Cont', 'Homo\nErr', 'Non-homo\nErr', 
                                  'Total\nErr', 'Anchor\nRate', 'Adj\nCont'],
                        cellLoc='center',
                        loc='center',
                        colColours=['#f0f8ff']*7)

single_table.auto_set_font_size(False)
single_table.set_fontsize(10)
single_table.scale(1.1, 1.8)

for j in range(7):
    single_table[(0, j)].set_facecolor('#4a86e8')
    single_table[(0, j)].set_text_props(color='white', weight='bold', fontsize=11)

ax1.set_title('Single-Resolution', fontsize=14, fontweight='bold', pad=20)

# 多分辨率表格 - 往上挪动
ax2.axis('tight')
ax2.axis('off')

multi_table_data = []
for _, row in df_multi.iterrows():
    multi_table_data.append([
        row['Method'],
        f"{row['Contiguity']:.4f}",
        f"{row['Inter_homo_error_rate']:.4f}",
        f"{row['Inter_nonhomo_error_rate']:.4f}",
        f"{row['total_error_rate']:.4f}",
        f"{row['Anchoring rate']:.4f}",
        f"{row['Adjusted contiguity']:.4f}"
    ])

multi_table = ax2.table(cellText=multi_table_data,
                       colLabels=['Method', 'Cont', 'Homo\nErr', 'Non-homo\nErr', 
                                 'Total\nErr', 'Anchor\nRate', 'Adj\nCont'],
                       cellLoc='center',
                       loc='center',
                       colColours=['#fff0f5']*7)

multi_table.auto_set_font_size(False)
multi_table.set_fontsize(10)
multi_table.scale(1.3, 1.5)  # 稍微调整行高

for j in range(7):
    multi_table[(0, j)].set_facecolor('#d9534f')
    multi_table[(0, j)].set_text_props(color='white', weight='bold', fontsize=11)

ax2.set_title('Multi-Resolution', fontsize=14, fontweight='bold', pad=0, y=0.8)  # 减少pad值让标题更靠近表格

# 调整子图间距，让多分辨率表格更靠上
plt.suptitle('C88 N50=2M - Compact Results Table', fontsize=16, fontweight='bold', y=0.98)
plt.tight_layout(rect=[0, 0, 1, 0.96], h_pad=0.5)  # 减小h_pad让两个表格更靠近

compact_filename = "C88_Compact_Results_Table_Fixed.png"
plt.savefig(compact_filename, dpi=300, bbox_inches='tight')
plt.show()

print(f"修正后的紧凑表格已保存为: {compact_filename}")

# 4. 创建一个更优的折线图，所有指标在一个图中
fig4, ax = plt.subplots(figsize=(20, 8))

# 准备数据
metrics_to_plot = [
    ('Adjusted contiguity', 'Adj Contiguity', 'blue', 'o-', 1.0),
    ('Anchoring rate', 'Anchoring Rate', 'green', 's-', 1.0),
    ('Contiguity', 'Contiguity', 'red', '^-', 100),  # Contiguity需要乘以100
    ('total_error_rate', 'Total Error', 'orange', 'd-', 1.0)
]

# 绘制所有指标
for metric_key, metric_name, color, marker, scale in metrics_to_plot:
    # 单分辨率数据
    single_values = df_single[metric_key].values * scale
    
    # 多分辨率数据
    multi_values = df_multi[metric_key].values * scale
    
    # 合并数据
    all_values = np.concatenate([single_values, multi_values])
    
    # 绘制折线
    ax.plot(x_positions, all_values, marker, color=color, linewidth=2, 
            markersize=8, label=metric_name, alpha=0.8)

# 添加垂直线区分单分辨率和多分辨率
ax.axvline(x=len(single_labels)-0.5, color='gray', linestyle='--', 
           alpha=0.7, linewidth=2, label='Single/Multi Boundary')

# 设置x轴
ax.set_xticks(x_positions)
ax.set_xticklabels(all_labels, rotation=45, ha='right', fontsize=10)

# 设置y轴
ax.set_ylabel('Score / Rate (%)', fontsize=12)

# 设置标题
ax.set_title('C88 N50=2M - All Performance Metrics Comparison', fontsize=14, fontweight='bold')

# 添加网格
ax.grid(True, alpha=0.3, linestyle='--')

# 添加图例
ax.legend(loc='upper right', fontsize=10, ncol=2)

# 添加区域标注
ax.text(len(single_labels)/2 - 1, ax.get_ylim()[1]*0.95, 'Single-Resolution', 
        ha='center', fontsize=12, fontweight='bold', bbox=dict(boxstyle='round,pad=0.3', facecolor='lightblue', alpha=0.5))
ax.text(len(single_labels) + len(multi_labels)/2 - 1, ax.get_ylim()[1]*0.95, 'Multi-Resolution', 
        ha='center', fontsize=12, fontweight='bold', bbox=dict(boxstyle='round,pad=0.3', facecolor='lightcoral', alpha=0.5))

plt.tight_layout()
combined_line_filename = "C88_All_Metrics_Combined_Line.png"
plt.savefig(combined_line_filename, dpi=300, bbox_inches='tight')
plt.show()

print(f"所有指标合并折线图已保存为: {combined_line_filename}")

print("\n" + "="*60)
print("生成的文件总结:")
print("="*60)
print("1. C88_Performance_Line_Charts.png - 4个子图的折线图")
print("2. C88_All_Results_Table_Fixed.png - 修正标红的完整表格")
print("3. C88_Compact_Results_Table_Fixed.png - 紧凑表格（多尺度表格上移）")
print("4. C88_All_Metrics_Combined_Line.png - 所有指标合并折线图")
print("="*60)