import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib.patches import Patch
from matplotlib.gridspec import GridSpec

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
    'HaphiC': {
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

# 指标名称
categories = [
    'Syntenic contigs\n(higher is better)',
    'Unanchored contigs\n(lower is better)', 
    'Newly anchored contigs\n(higher is better)',
    'Translocation contigs\n(lower is better)',
    'Relocation contigs\n(lower is better)',
    'Inversion contigs\n(lower is better)',
    'Inversion and relocation\ncontigs (lower is better)'
]

# 工具名称
tools = ['ALLHiC', 'HaphiC', 'ConHiC (cluster)', 'ConHiC (cluster+sort)']

# 数据集名称
datasets = ['C88(2M)', 'C88(0.8M)', 'SS', 'XJDY', 'ZM-4', 'cs', 'ap']

# 颜色设置 - 使用学术常用的成对配色
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']  # 蓝、橙、绿、红

# 创建保存目录
output_dir = 'bar_charts'
os.makedirs(output_dir, exist_ok=True)

# 设置全局样式 - 使用系统通用字体
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans', 'Bitstream Vera Sans', 'sans-serif']
plt.rcParams['font.size'] = 11  # 修改: 从9改为11，确保大于10
plt.rcParams['axes.linewidth'] = 0.8
# EPS格式相关设置
plt.rcParams['ps.useafm'] = True
plt.rcParams['pdf.use14corefonts'] = True
plt.rcParams['text.usetex'] = False  # 不使用LaTeX以避免兼容性问题

print("正在生成柱状图（纵坐标从最小值的80%开始，所有字体≥10）...")

# 创建7个子图，布局为2行+2行+3行
fig = plt.figure(figsize=(20, 16))  # 增加图形尺寸以适应更大的字体

# 定义子图网格
# 第一行2个，第二行2个，第三行3个
ax1 = plt.subplot(3, 2, 1)
ax2 = plt.subplot(3, 2, 2)
ax3 = plt.subplot(3, 2, 3)
ax4 = plt.subplot(3, 2, 4)
# 第三行使用特殊的gridspec来创建3列
gs = GridSpec(3, 3, figure=fig)
ax5 = fig.add_subplot(gs[2, 0])
ax6 = fig.add_subplot(gs[2, 1])
ax7 = fig.add_subplot(gs[2, 2])

axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7]

# 柱状图宽度和位置设置
bar_width = 0.2
x = np.arange(len(datasets))

# 为每个指标创建柱状图
for metric_idx, (category, ax) in enumerate(zip(categories, axes)):
    
    # 为每个工具收集数据
    tool_data = []
    has_data = []  # 记录每个工具在该指标上是否有数据
    
    for tool_idx, tool in enumerate(tools):
        tool_metric_data = []
        tool_has_data = []
        for dataset in datasets:
            dataset_data = data[tool][dataset]
            if dataset_data is not None:
                value = dataset_data[metric_idx]
                tool_metric_data.append(value)
                tool_has_data.append(True)
            else:
                tool_metric_data.append(0)  # 填充0，但不会显示
                tool_has_data.append(False)
        tool_data.append(tool_metric_data)
        has_data.append(tool_has_data)
    
    # 绘制柱状图
    for tool_idx in range(len(tools)):
        # 计算每个柱的位置
        bars = ax.bar(x + tool_idx * bar_width, tool_data[tool_idx], 
                      bar_width, label=tools[tool_idx] if metric_idx == 0 else "",
                      color=colors[tool_idx], alpha=0.8, edgecolor='black', linewidth=0.5)
        
        # 在柱子上添加数值标签
        for i, (value, has_val) in enumerate(zip(tool_data[tool_idx], has_data[tool_idx])):
            if has_val and value > 0:  # 只显示有数据的柱子
                height = value
                # 修改: fontsize从6改为10
                ax.text(i + tool_idx * bar_width, height + 0.5, f'{value:.2f}', 
                        ha='center', va='bottom', fontsize=10, rotation=90)
    
    # 设置子图标题和标签
    ax.set_title(category, fontsize=12, pad=15)  # 修改: 从11改为12
    ax.set_ylabel('Percentage (%)', fontsize=11)  # 修改: 从9改为11
    
    # 修改: x轴标签字体大小从8改为10
    ax.set_xticks(x + bar_width * 1.5)  # 将刻度放在四个柱子的中间
    
    # 旋转x轴标签以避免重叠
    if metric_idx in [4, 5, 6]:  # 最后一行子图
        ax.set_xticklabels(datasets, fontsize=10, rotation=45, ha='right')  # 修改: 从8改为10
    else:
        ax.set_xticklabels(datasets, fontsize=10)  # 修改: 从8改为10
    
    # 计算纵坐标范围 - 从最小值的80%开始
    # 收集所有有效数据
    all_valid_values = []
    for tool_data_i in tool_data:
        for val in tool_data_i:
            if val > 0:  # 只考虑有效值（大于0）
                all_valid_values.append(val)
    
    if all_valid_values:
        min_val = min(all_valid_values)
        max_val = max(all_valid_values)
        
        # 设置y轴从最小值的80%开始
        y_min = min_val * 0.8
        y_max = max_val * 1.2  # 修改: 从1.15增加到1.2，为更大的数值标签留出更多空间
        
        # 确保y_min不小于0
        y_min = max(0, y_min)
        
        ax.set_ylim(y_min, y_max)
        
        # 调整y轴刻度，使其从y_min开始
        # 计算合适的刻度间隔
        range_val = y_max - y_min
        if range_val > 20:
            step = 5
        elif range_val > 10:
            step = 2
        elif range_val > 5:
            step = 1
        elif range_val > 2:
            step = 0.5
        else:
            step = 0.2
        
        # 生成从y_min开始的刻度
        ticks = np.arange(np.ceil(y_min/step)*step, y_max + step, step)
        ax.set_yticks(ticks)
        
        # 格式化y轴刻度标签，保留1-2位小数
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x:.1f}' if x >= 10 else f'{x:.2f}'))
    else:
        # 如果没有有效数据，使用默认范围
        ax.set_ylim(0, 10)
    
    # 添加网格线
    ax.yaxis.grid(True, linestyle='--', alpha=0.3, linewidth=0.5)
    ax.set_axisbelow(True)
    
    # 修改: y轴刻度字体大小从8改为10
    ax.tick_params(axis='y', labelsize=10)
    
    # 为需要显示"lower is better"的指标添加特殊标记
    # 修改: fontsize从7改为10
    if metric_idx in [1, 3, 4, 5, 6]:  # 这些指标是lower is better
        # 在标题或角落添加标记
        ax.text(0.02, 0.95, '↓ lower is better', transform=ax.transAxes, 
                fontsize=10, verticalalignment='top', color='#555555', style='italic')
    else:
        ax.text(0.02, 0.95, '↑ higher is better', transform=ax.transAxes, 
                fontsize=10, verticalalignment='top', color='#555555', style='italic')

# 添加图例
handles = [Patch(color=colors[i], label=tools[i]) for i in range(len(tools))]
fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, 0.985), 
           ncol=4, fontsize=11, frameon=False)  # 修改: 从10改为11

# 调整布局
plt.tight_layout()
plt.subplots_adjust(top=0.93, hspace=0.4, wspace=0.3, bottom=0.1)  # 调整间距以适应更大的字体

# 保存图片 - PNG格式
output_path = os.path.join(output_dir, 'performance_bar_charts_7metrics_zoomed_font10plus.png')
plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
print(f"柱状图（缩放版，字体≥10）已保存至: {output_path}")

# 保存PDF格式
output_path_pdf = os.path.join(output_dir, 'performance_bar_charts_7metrics_zoomed_font10plus.pdf')
plt.savefig(output_path_pdf, bbox_inches='tight', facecolor='white')
print(f"PDF版本已保存至: {output_path_pdf}")

# 保存EPS格式
output_path_eps = os.path.join(output_dir, 'performance_bar_charts_7metrics_zoomed_font10plus.eps')
plt.savefig(output_path_eps, format='eps', bbox_inches='tight', facecolor='white', dpi=300)
print(f"EPS版本已保存至: {output_path_eps}")

plt.show()

# 打印字体大小统计
print(f"\n{'='*70}")
print("字体大小统计（所有字体≥10）:")
print(f"{'='*70}")
print(f"全局字体: {plt.rcParams['font.size']}")
print(f"子图标题: 12")
print(f"y轴标签: 11")
print(f"x轴标签: 10")
print(f"y轴刻度: 10")
print(f"数值标签: 10")
print(f"图例: 11")
print(f"标记文字: 10")
print(f"{'='*70}")

# 打印每个指标的纵坐标范围信息
print(f"\n{'='*70}")
print("各指标的纵坐标范围（从最小值的80%开始）:")
print(f"{'='*70}")

for metric_idx, category in enumerate(categories):
    # 收集所有有效数据
    all_valid_values = []
    for tool_idx, tool in enumerate(tools):
        for dataset in datasets:
            dataset_data = data[tool][dataset]
            if dataset_data is not None:
                value = dataset_data[metric_idx]
                if value > 0:
                    all_valid_values.append(value)
    
    if all_valid_values:
        min_val = min(all_valid_values)
        max_val = max(all_valid_values)
        y_min = max(0, min_val * 0.8)
        y_max = max_val * 1.2  # 修改为1.2以匹配上面的调整
        
        print(f"\n{category}:")
        print(f"  数据范围: [{min_val:.3f}, {max_val:.3f}]")
        print(f"  纵坐标范围: [{y_min:.3f}, {y_max:.3f}]")
        print(f"  缩放因子: {(max_val - y_min)/(max_val - min_val):.2f}x")

print(f"\n{'='*70}")
print(f"所有图表已成功保存至目录: {output_dir}")
print(f"输出文件:")
print(f"  - {output_path}")
print(f"  - {output_path_pdf}")
print(f"  - {output_path_eps}")
print(f"{'='*70}")