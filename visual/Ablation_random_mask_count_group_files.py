import os
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

def collect_data_for_resolution(base_path, resolution_name):
    """
    收集指定路径和分辨率的数据，包括group files和main groups
    """
    base_dir = Path(base_path)
    if not base_dir.exists():
        print(f"警告: 路径不存在 {base_dir}")
        return None, None, None, None
    
    # 收集数据
    table_data = []  # group files数据
    main_groups_data = []  # main groups数据
    all_counts = []  # 所有group files计数
    all_main_groups = []  # 所有main groups计数
    
    for random_num in range(10):
        row_data = []
        main_row_data = []
        for inflation in [1.1, 1.2, 1.3]:
            inflation_dir = base_dir / str(random_num) / f"inflation_{inflation}"
            group_count = 0
            main_groups = 0
            
            if inflation_dir.exists():
                # 统计group文件数量
                group_files = list(inflation_dir.glob("group*"))
                group_count = len(group_files)
                
                # 尝试读取main groups数据
                main_groups_file = inflation_dir / "main_groups.txt"
                if main_groups_file.exists():
                    try:
                        with open(main_groups_file, 'r') as f:
                            content = f.read().strip()
                            if content:
                                main_groups = int(content)
                    except:
                        main_groups = 0
                else:
                    main_groups = 0
            
            row_data.append(group_count)
            main_row_data.append(main_groups)
            all_counts.append(group_count)
            all_main_groups.append(main_groups)
        
        table_data.append(row_data)
        main_groups_data.append(main_row_data)
    
    return table_data, main_groups_data, all_counts, all_main_groups

def create_combined_merged_table():
    """
    创建包含三个分辨率的合并单元格表格图
    """
    # 1. 收集三个分辨率的数据
    print("开始收集数据...")
    
    # mean/10 数据  
    path_mean10 = "/home/chenzh/HapHiC/data/C88/Random_mask/test28/01.cluster/108582"
    data_mean10, main_mean10, counts_mean10, main_counts_mean10 = collect_data_for_resolution(path_mean10, "mean/10")
    
    # mean/8 数据  
    path_mean8 = "/home/chenzh/HapHiC/data/C88/Random_mask/test28/01.cluster/135727"
    data_mean8, main_mean8, counts_mean8, main_counts_mean8 = collect_data_for_resolution(path_mean8, "mean/8")
    
    # mean/6 数据  
    path_mean6 = "/home/chenzh/HapHiC/data/C88/Random_mask/test28/01.cluster/180970"
    data_mean6, main_mean6, counts_mean6, main_counts_mean6 = collect_data_for_resolution(path_mean6, "mean/6")
    
    # 检查数据是否完整
    resolutions_data = {
        'mean/10': (data_mean10, main_mean10, counts_mean10, main_counts_mean10, path_mean10),
        'mean/8': (data_mean8, main_mean8, counts_mean8, main_counts_mean8, path_mean8),
        'mean/6': (data_mean6, main_mean6, counts_mean6, main_counts_mean6, path_mean6)
    }
    
    for res_name, (data, main_data, counts, main_counts, path) in resolutions_data.items():
        if data is None:
            print(f"错误: 无法读取{res_name}数据，请检查路径: {path}")
            return None
    
    print(f"✓ mean/10: 找到 {len(counts_mean10)} 个数据点, {len(main_counts_mean10)} 个main groups")
    print(f"✓ mean/8: 找到 {len(counts_mean8)} 个数据点, {len(main_counts_mean8)} 个main groups")
    print(f"✓ mean/6: 找到 {len(counts_mean6)} 个数据点, {len(main_counts_mean6)} 个main groups")
    
    # 2. 使用提供的main groups数据覆盖默认值
    # mean/10 数据
    main_data_mean10 = [
        [46, 824, 1058],  # random 0
        [45, 1, 1056],    # random 1
        [46, 829, 1056],  # random 2
        [46, 1, 1056],    # random 3
        [1, 816, 1056],   # random 4
        [47, 814, 1057],  # random 5
        [46, 824, 1058],  # random 6
        [46, 1, 1058],    # random 7
        [45, 835, 1059],  # random 8
        [44, 809, 1055]   # random 9
    ]
    
    # mean/8 数据
    main_data_mean8 = [
        [1, 730, 1054],   # random 0
        [2, 748, 1054],   # random 1
        [44, 727, 1055],  # random 2
        [1, 725, 1051],   # random 3
        [47, 1, 1054],    # random 4
        [44, 737, 1054],  # random 5
        [44, 725, 1052],  # random 6
        [44, 709, 1054],  # random 7
        [44, 736, 1049],  # random 8
        [44, 730, 1053]   # random 9
    ]
    
    # mean/6 数据
    main_data_mean6 = [
        [1, 620, 1044],   # random 0
        [1, 592, 1050],   # random 1
        [2, 629, 1043],   # random 2
        [2, 606, 1050],   # random 3
        [3, 619, 1047],   # random 4
        [3, 603, 1044],   # random 5
        [1, 612, 1045],   # random 6
        [1, 1, 1046],     # random 7
        [1, 613, 1047],   # random 8
        [1, 613, 1045]    # random 9
    ]
    
    # 更新main groups数据
    main_mean10 = main_data_mean10
    main_mean8 = main_data_mean8
    main_mean6 = main_data_mean6
    
    # 更新main groups计数
    main_counts_mean10 = [item for sublist in main_mean10 for item in sublist]
    main_counts_mean8 = [item for sublist in main_mean8 for item in sublist]
    main_counts_mean6 = [item for sublist in main_mean6 for item in sublist]
    
    # 3. 创建图形和坐标轴 - 增大图形尺寸以容纳更大字号
    fig, ax = plt.subplots(figsize=(18, 14))  # 增大图形尺寸
    ax.axis('off')
    
    # 4. 表格参数设置
    # 每个分辨率有30行数据 + 1行表头
    n_data_rows_per_res = 30  # 10个random × 3个inflation
    n_resolutions = 3  # mean/10, mean/8, mean/6
    n_total_rows = n_data_rows_per_res * n_resolutions + 1  # 加上表头行
    n_cols = 5  # 5列: Resolution, Random number, Inflation, Group number, Main groups
    
    # 调整单元格尺寸以容纳更大字号
    cell_width = 0.16  # 稍微增加宽度
    cell_height = 0.032  # 增加高度以容纳更大字号
    table_left = 0.10
    table_bottom = 0.06
    
    # 5. 绘制表头 - 增大表头字号
    headers = ['Resolution', 'Random number', 'Inflation', 'Group number', 'Main groups']
    header_colors = ['#4CAF50'] * n_cols
    
    # 表头位置：最上面一行
    header_y = table_bottom + (n_total_rows - 1) * cell_height
    
    for col in range(n_cols):
        x = table_left + col * cell_width
        
        # 绘制表头单元格
        rect = Rectangle((x, header_y), cell_width, cell_height,
                         facecolor=header_colors[col], edgecolor='black', linewidth=1.5)  # 增加边框线宽
        ax.add_patch(rect)
        
        # 添加表头文本 - 增大字号
        ax.text(x + cell_width/2, header_y + cell_height/2, headers[col],
                ha='center', va='center', color='white', fontweight='bold', fontsize=14)  # 从9增大到14
    
    # 6. 绘制第一列：Resolution合并单元格
    resolution_height = n_data_rows_per_res * cell_height
    
    # 定义分辨率顺序和颜色
    resolutions = [
        ('mean/10', '#E8F8E8'),  # 浅绿色
        ('mean/8', '#FFF0F0'),   # 浅红色
        ('mean/6', '#FFF8E1')    # 浅黄色
    ]
    
    # 为每个分辨率绘制合并单元格
    for i, (res_name, color) in enumerate(resolutions):
        y_pos = table_bottom + i * resolution_height
        
        # 绘制分辨率单元格
        rect_res = Rectangle((table_left, y_pos), 
                            cell_width, resolution_height,
                            facecolor=color, edgecolor='black', linewidth=1.5)
        ax.add_patch(rect_res)
        
        # 计算中心位置并添加文本 - 增大字号
        center_y = y_pos + resolution_height / 2
        ax.text(table_left + cell_width/2, center_y, res_name,
                ha='center', va='center', fontsize=13, fontweight='bold', rotation=0)  # 从9增大到13
    
    # 7. 绘制第二列：Random number（每3行合并）
    for res_idx in range(n_resolutions):
        y_base = table_bottom + res_idx * resolution_height
        
        for random_num in range(10):
            block_start_row = random_num * 3
            y_pos = y_base + (n_data_rows_per_res - block_start_row - 3) * cell_height
            
            rect_random = Rectangle((table_left + cell_width, y_pos), 
                                   cell_width, cell_height * 3,
                                   facecolor='#F5F5F5', edgecolor='black', linewidth=1.5)
            ax.add_patch(rect_random)
            
            random_center_y = y_pos + (cell_height * 3) / 2
            ax.text(table_left + cell_width + cell_width/2, random_center_y, str(random_num),
                    ha='center', va='center', fontsize=12, fontweight='bold')  # 从7增大到12
    
    # 8. 计算最大值用于颜色标准化
    max_values = {
        'mean/10': max(counts_mean10) if counts_mean10 else 1,
        'mean/8': max(counts_mean8) if counts_mean8 else 1,
        'mean/6': max(counts_mean6) if counts_mean6 else 1
    }
    
    max_main_values = {
        'mean/10': max(main_counts_mean10) if main_counts_mean10 else 1,
        'mean/8': max(main_counts_mean8) if main_counts_mean8 else 1,
        'mean/6': max(main_counts_mean6) if main_counts_mean6 else 1
    }
    
    # 定义每个分辨率的数据
    resolution_info = [
        ('mean/10', data_mean10, main_mean10, max_values['mean/10'], max_main_values['mean/10'], (0.8, 1.0, 0.9)),
        ('mean/8', data_mean8, main_mean8, max_values['mean/8'], max_main_values['mean/8'], (1.0, 0.9, 0.9)),
        ('mean/6', data_mean6, main_mean6, max_values['mean/6'], max_main_values['mean/6'], (1.0, 0.95, 0.8))
    ]
    
    # 9. 绘制第三、四、五列数据
    for res_idx, (res_name, data, main_data, max_val, max_main_val, base_color) in enumerate(resolution_info):
        y_base = table_bottom + res_idx * resolution_height
        data_row_idx = 0
        
        for random_num in range(10):
            for infl_idx, inflation in enumerate([1.1, 1.2, 1.3]):
                y_pos = y_base + (n_data_rows_per_res - data_row_idx - 1) * cell_height
                
                # 第三列: Inflation - 增大字号
                rect_infl = Rectangle((table_left + 2 * cell_width, y_pos), 
                                     cell_width, cell_height,
                                     facecolor='#FFFFFF', edgecolor='black', linewidth=1.5)
                ax.add_patch(rect_infl)
                ax.text(table_left + 2 * cell_width + cell_width/2, y_pos + cell_height/2, 
                       str(inflation), ha='center', va='center', fontsize=11)  # 从7增大到11
                
                # 第四列: Group number - 增大字号
                count_value = data[random_num][infl_idx]
                
                if max_val > 0:
                    color_intensity = 0.3 + count_value / max_val * 0.6
                    color_intensity = min(0.9, color_intensity)
                else:
                    color_intensity = 0.3
                
                if res_name == 'mean/10':
                    cell_color = (base_color[0] - color_intensity * 0.1, base_color[1], base_color[2] - color_intensity * 0.1)
                elif res_name == 'mean/8':
                    cell_color = (base_color[0], base_color[1] - color_intensity * 0.1, base_color[2] - color_intensity * 0.1)
                else:
                    cell_color = (base_color[0], base_color[1] - color_intensity * 0.1, base_color[2] - color_intensity * 0.2)
                
                rect_count = Rectangle((table_left + 3 * cell_width, y_pos), 
                                      cell_width, cell_height,
                                      facecolor=cell_color, edgecolor='black', linewidth=1.5)
                ax.add_patch(rect_count)
                ax.text(table_left + 3 * cell_width + cell_width/2, y_pos + cell_height/2, 
                       str(count_value), ha='center', va='center', fontsize=11, fontweight='bold')  # 从7增大到11
                
                # 第五列: Main groups - 增大字号
                main_value = main_data[random_num][infl_idx]
                
                if main_value > 1000:
                    main_cell_color = (0.2, 0.4, 0.8)
                elif main_value > 500:
                    main_cell_color = (0.4, 0.6, 0.9)
                elif main_value > 100:
                    main_cell_color = (0.6, 0.8, 1.0)
                elif main_value > 10:
                    main_cell_color = (0.8, 0.9, 1.0)
                else:
                    main_cell_color = (1.0, 1.0, 1.0)
                
                rect_main = Rectangle((table_left + 4 * cell_width, y_pos), 
                                     cell_width, cell_height,
                                     facecolor=main_cell_color, edgecolor='black', linewidth=1.5)
                ax.add_patch(rect_main)
                ax.text(table_left + 4 * cell_width + cell_width/2, y_pos + cell_height/2, 
                       str(main_value), ha='center', va='center', fontsize=11, fontweight='bold')  # 从7增大到11
                
                data_row_idx += 1
    
    # 10. 设置图形边界
    table_width = n_cols * cell_width
    table_height = n_total_rows * cell_height
    ax.set_xlim(table_left - 0.05, table_left + table_width + 0.05)
    ax.set_ylim(table_bottom - 0.05, table_bottom + table_height + 0.05)
    
    # 11. 添加标题和统计信息 - 增大标题字号
    total_counts = {
        'mean/10': sum(counts_mean10),
        'mean/8': sum(counts_mean8),
        'mean/6': sum(counts_mean6)
    }
    
    total_main_counts = {
        'mean/10': sum(main_counts_mean10),
        'mean/8': sum(main_counts_mean8),
        'mean/6': sum(main_counts_mean6)
    }
    
    avg_counts = {
        'mean/10': total_counts['mean/10'] / len(counts_mean10) if counts_mean10 else 0,
        'mean/8': total_counts['mean/8'] / len(counts_mean8) if counts_mean8 else 0,
        'mean/6': total_counts['mean/6'] / len(counts_mean6) if counts_mean6 else 0
    }
    
    avg_main_counts = {
        'mean/10': total_main_counts['mean/10'] / len(main_counts_mean10) if main_counts_mean10 else 0,
        'mean/8': total_main_counts['mean/8'] / len(main_counts_mean8) if main_counts_mean8 else 0,
        'mean/6': total_main_counts['mean/6'] / len(main_counts_mean6) if main_counts_mean6 else 0
    }
    
    title_text = 'Group File Counts and Main Groups - Three Resolutions\n'
    title_text += f"mean/10: {total_counts['mean/10']} groups (avg: {avg_counts['mean/10']:.1f}), "
    title_text += f"{total_main_counts['mean/10']} main groups (avg: {avg_main_counts['mean/10']:.1f})\n"
    title_text += f"mean/8: {total_counts['mean/8']} groups (avg: {avg_counts['mean/8']:.1f}), "
    title_text += f"{total_main_counts['mean/8']} main groups (avg: {avg_main_counts['mean/8']:.1f})\n"
    title_text += f"mean/6: {total_counts['mean/6']} groups (avg: {avg_counts['mean/6']:.1f}), "
    title_text += f"{total_main_counts['mean/6']} main groups (avg: {avg_main_counts['mean/6']:.1f})"
    
    plt.title(title_text, fontsize=13, fontweight='bold', pad=30, loc='center')  # 从10增大到13
    
    # 12. 添加分辨率之间的分隔线
    for i in range(1, n_resolutions):
        separator_y = table_bottom + i * resolution_height
        ax.axhline(y=separator_y, xmin=table_left, xmax=table_left + table_width, 
                  color='red', linestyle='--', linewidth=1.5, alpha=0.5)
    
    # 13. 保存和显示
    output_file = "group_counts_three_resolutions_with_main.png"
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\n✓ 三分辨率表格图已保存到: {output_file}")
    print(f"✓ 表格尺寸: {n_total_rows}行 × {n_cols}列")
    print(f"✓ 包含列: Resolution, Random number, Inflation, Group number, Main groups")
    
    # 显示图形
    plt.show()
    
    return output_file, data_mean10, main_mean10, data_mean8, main_mean8, data_mean6, main_mean6

def print_data_summary(data_mean10, main_mean10, data_mean8, main_mean8, data_mean6, main_mean6):
    """打印数据摘要"""
    print("\n" + "="*80)
    print("数据摘要 - 三个分辨率 (包含Main groups)")
    print("="*80)
    
    resolutions = [
        ('mean/10', data_mean10, main_mean10),
        ('mean/8', data_mean8, main_mean8),
        ('mean/6', data_mean6, main_mean6)
    ]
    
    total_all = 0
    total_main_all = 0
    
    for res_name, data, main_data in resolutions:
        print(f"\n{res_name}:")
        print("-"*50)
        total = 0
        total_main = 0
        
        for i in range(10):
            row = data[i]
            main_row = main_data[i]
            row_sum = sum(row)
            main_row_sum = sum(main_row)
            total += row_sum
            total_main += main_row_sum
            print(f"Random {i}: Groups [1.1:{row[0]:3d}, 1.2:{row[1]:3d}, 1.3:{row[2]:3d}] Total: {row_sum:4d} | "
                  f"Main [1.1:{main_row[0]:4d}, 1.2:{main_row[1]:4d}, 1.3:{main_row[2]:4d}] Total: {main_row_sum:5d}")
        
        print(f"总计: {total} 个group文件, {total_main} 个main groups")
        avg = total / 30 if total > 0 else 0
        avg_main = total_main / 30 if total_main > 0 else 0
        print(f"平均每目录: {avg:.1f} 个groups, {avg_main:.1f} 个main groups")
        
        total_all += total
        total_main_all += total_main
    
    print("\n" + "="*80)
    print(f"所有分辨率总计: {total_all} 个group文件, {total_main_all} 个main groups")
    print(f"总平均每目录: {total_all/90:.1f} 个groups, {total_main_all/90:.1f} 个main groups")

def main():
    print("开始创建三分辨率合并单元格表格图(包含Main groups)...")
    print("="*80)
    
    # 创建包含三个分辨率的表格图
    result = create_combined_merged_table()
    
    if result:
        output_file, data_mean10, main_mean10, data_mean8, main_mean8, data_mean6, main_mean6 = result
        # 打印数据摘要
        print_data_summary(data_mean10, main_mean10, data_mean8, main_mean8, data_mean6, main_mean6)
        
        print(f"\n✓ 表格生成完成!")
        print(f"✓ 输出文件: {output_file}")

if __name__ == "__main__":
    main()