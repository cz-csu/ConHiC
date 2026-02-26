import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.colors as mcolors

# 获取Paired配色方案
Paired_3 = plt.cm.Paired
colors = Paired_3.colors      # 成对配色 (学术常用)

def create_table1():
    """创建第一个表格：C88 N50=2M 分辨率对比（使用缩写）"""
    fig, ax = plt.subplots(figsize=(10, 4.5))
    ax.axis('off')
    
    # 表格数据 - 使用缩写的表头
    headers = ['Resolution', 'Contiguity', 'Homo_err', 
               'Nonhomo_err', 'Total_err', 
               'Anchor_rate', 'Adj_contig', 'Clusters']
    
    data = [
        ['mean/6', '-', '-', '-', '-', '-', '-', '-'],
        ['mean/8', '0.8552', '7.8676', '5.0622', '12.9298', '87.5251', '74.8515', '48'],
        ['mean/10', '0.7736', '9.8756', '10.6296', '20.5052', '82.3276', '63.6886', '48']
    ]
    
    # 表格参数
    n_rows = len(data) + 1  # 数据行 + 表头
    n_cols = len(headers)
    cell_height = 0.15
    cell_width = 0.105
    table_left = 0.07
    table_bottom = 0.15
    
    # 绘制表头 - 使用深色
    header_y = table_bottom + (n_rows - 1) * cell_height
    for col in range(n_cols):
        x = table_left + col * cell_width
        rect = Rectangle((x, header_y), cell_width, cell_height,
                         facecolor=colors[5], edgecolor='black', linewidth=1.2)  # 使用深蓝色
        ax.add_patch(rect)
        ax.text(x + cell_width/2, header_y + cell_height/2, headers[col],
                ha='center', va='center', color='white', fontweight='bold', fontsize=10)
    
    # 绘制数据行
    for row_idx, row_data in enumerate(data):
        y_pos = table_bottom + (n_rows - row_idx - 2) * cell_height
        
        # 使用Paired配色交替行颜色
        if row_idx % 2 == 0:
            row_color = colors[0]  # 浅蓝色
        else:
            row_color = colors[1]  # 浅橙色
        
        for col in range(n_cols):
            x = table_left + col * cell_width
            rect = Rectangle((x, y_pos), cell_width, cell_height,
                            facecolor=row_color, edgecolor='black', linewidth=1.0,
                            alpha=0.7)  # 添加透明度使颜色更柔和
            ax.add_patch(rect)
            
            cell_value = row_data[col]
            ax.text(x + cell_width/2, y_pos + cell_height/2, cell_value,
                    ha='center', va='center', color='black', 
                    fontweight='normal', fontsize=10)
    
    # 设置图形边界
    table_width = n_cols * cell_width
    table_height = n_rows * cell_height
    ax.set_xlim(table_left - 0.02, table_left + table_width + 0.02)
    ax.set_ylim(table_bottom - 0.02, table_bottom + table_height + 0.05)
    
    # 添加标题
    plt.title('Table 1: C88 N50=2M - Resolution Comparison', 
              fontsize=13, fontweight='bold', pad=15, color='black')
    
    plt.tight_layout()
    plt.savefig('table1_resolution_comparison_paired.png', dpi=300, bbox_inches='tight')
    print("✓ Table 1 saved as 'table1_resolution_comparison_paired.png'")
    plt.show()

def create_table2():
    """创建第二个表格：Random mask效果对比（使用缩写）"""
    fig, ax = plt.subplots(figsize=(11, 4.5))
    ax.axis('off')
    
    # 表格数据 - 使用缩写的表头
    headers = ['Method', 'Contiguity', 'Homo_err', 
               'Nonhomo_err', 'Total_err', 
               'Anchor_rate', 'Adj_contig']
    
    data = [
        ['Random mask', '0.8732', '5.4975', '4.9570', '10.4545', '89.9603', '78.5582'],
        ['w/o Random mask', '0.8559', '9.3436', '8.0841', '17.4277', '89.5845', '76.6754'],
        ['w/o Random mask\n(no mean/6)', '0.8663', '6.0501', '4.8885', '10.9386', '89.6536', '77.6669']
    ]
    
    # 表格参数
    n_rows = len(data) + 1
    n_cols = len(headers)
    cell_height = 0.16
    cell_width = 0.125
    table_left = 0.07
    table_bottom = 0.15
    
    # 绘制表头 - 使用深色
    header_y = table_bottom + (n_rows - 1) * cell_height
    for col in range(n_cols):
        x = table_left + col * cell_width
        rect = Rectangle((x, header_y), cell_width, cell_height,
                         facecolor=colors[5], edgecolor='black', linewidth=1.2)  # 使用深蓝色
        ax.add_patch(rect)
        ax.text(x + cell_width/2, header_y + cell_height/2, headers[col],
                ha='center', va='center', color='white', fontweight='bold', fontsize=11)
    
    # 绘制数据行
    for row_idx, row_data in enumerate(data):
        y_pos = table_bottom + (n_rows - row_idx - 2) * cell_height
        
        # 使用Paired配色交替行颜色
        if row_idx % 2 == 0:
            row_color = colors[2]  # 浅绿色
        else:
            row_color = colors[3]  # 浅红色
        
        for col in range(n_cols):
            x = table_left + col * cell_width
            rect = Rectangle((x, y_pos), cell_width, cell_height,
                            facecolor=row_color, edgecolor='black', linewidth=1.0,
                            alpha=0.7)  # 添加透明度使颜色更柔和
            ax.add_patch(rect)
            
            cell_value = row_data[col]
            
            # 对于第一列（方法名称），使用加粗字体
            if col == 0:
                font_weight = 'bold'
                fontsize = 9
            else:
                font_weight = 'normal'
                fontsize = 10
            
            ax.text(x + cell_width/2, y_pos + cell_height/2, cell_value,
                    ha='center', va='center', color='black', 
                    fontweight=font_weight, fontsize=fontsize)
    
    # 设置图形边界
    table_width = n_cols * cell_width
    table_height = n_rows * cell_height
    ax.set_xlim(table_left - 0.02, table_left + table_width + 0.02)
    ax.set_ylim(table_bottom - 0.02, table_bottom + table_height + 0.05)
    
    # 添加标题
    plt.title('Table 2: C88(N50=2M) - Random Mask Effect Comparison', 
              fontsize=13, fontweight='bold', pad=15, color='black')
    
    plt.tight_layout()
    plt.savefig('table2_random_mask_comparison_paired.png', dpi=300, bbox_inches='tight')
    print("✓ Table 2 saved as 'table2_random_mask_comparison_paired.png'")
    plt.show()

def create_table3():
    """创建第三个表格：Contig类型统计（使用缩写）"""
    fig, ax = plt.subplots(figsize=(13, 4.5))
    ax.axis('off')
    
    # 表格数据 - 使用缩写的表头
    headers = ['Method', 'Syn_contigs', 'Unanchored', 
               'New_anchor', 'Transloc', 
               'Reloc', 'Inversion', 
               'Inv+Reloc']
    
    data = [
        ['Random mask', '57.2619', '8.8432', '1.1963', '9.4040', 
         '9.9290', '8.7165', '4.6479'],
        ['w/o Random mask', '51.7292', '9.2191', '1.1963', '17.6672', 
         '8.5821', '8.0442', '3.5615'],
        ['w/o Random mask\n(no mean/6)', '56.6051', '9.1500', '1.1963', 
         '9.8070', '9.6878', '9.0775', '4.4760']
    ]
    
    # 表格参数
    n_rows = len(data) + 1
    n_cols = len(headers)
    cell_height = 0.16
    cell_width = 0.11
    table_left = 0.06
    table_bottom = 0.15
    
    # 绘制表头 - 使用深色
    header_y = table_bottom + (n_rows - 1) * cell_height
    for col in range(n_cols):
        x = table_left + col * cell_width
        rect = Rectangle((x, header_y), cell_width, cell_height,
                         facecolor=colors[7], edgecolor='black', linewidth=1.2)  # 使用深紫色
        ax.add_patch(rect)
        
        fontsize = 9 if len(headers[col]) > 8 else 10
        ax.text(x + cell_width/2, header_y + cell_height/2, headers[col],
                ha='center', va='center', color='white', 
                fontweight='bold', fontsize=fontsize)
    
    # 绘制数据行
    for row_idx, row_data in enumerate(data):
        y_pos = table_bottom + (n_rows - row_idx - 2) * cell_height
        
        # 使用Paired配色交替行颜色
        if row_idx % 2 == 0:
            row_color = colors[4]  # 浅紫色
        else:
            row_color = colors[5]  # 浅蓝色
        
        for col in range(n_cols):
            x = table_left + col * cell_width
            rect = Rectangle((x, y_pos), cell_width, cell_height,
                            facecolor=row_color, edgecolor='black', linewidth=1.0,
                            alpha=0.7)  # 添加透明度使颜色更柔和
            ax.add_patch(rect)
            
            cell_value = row_data[col]
            
            # 对于第一列（方法名称），使用加粗字体
            if col == 0:
                font_weight = 'bold'
                fontsize = 9
            else:
                font_weight = 'normal'
                fontsize = 9.5
            
            ax.text(x + cell_width/2, y_pos + cell_height/2, cell_value,
                    ha='center', va='center', color='black', 
                    fontweight=font_weight, fontsize=fontsize)
    
    # 设置图形边界
    table_width = n_cols * cell_width
    table_height = n_rows * cell_height
    ax.set_xlim(table_left - 0.02, table_left + table_width + 0.02)
    ax.set_ylim(table_bottom - 0.02, table_bottom + table_height + 0.05)
    
    # 添加标题
    plt.title('Table 3: C88(N50=2M) - Contig Type Statistics (%)', 
              fontsize=13, fontweight='bold', pad=15, color='black')
    
    plt.tight_layout()
    plt.savefig('table3_contig_statistics_paired.png', dpi=300, bbox_inches='tight')
    print("✓ Table 3 saved as 'table3_contig_statistics_paired.png'")
    plt.show()

def create_combined_tables():
    """创建包含两个表格的汇总图（使用Paired配色）- 去掉子图(c)"""
    fig, axes = plt.subplots(2, 1, figsize=(13, 10))  # 改为2行1列，调整figsize高度
    #fig.suptitle('C88 N50=2M - Comprehensive Analysis', 
    #             fontsize=15, fontweight='bold', y=0.98, color='black')
    
    # 隐藏所有坐标轴
    for ax in axes:
        ax.axis('off')
    
    # 表格1数据（进一步缩写）
    headers1 = ['Resolution', 'Contiguity', 'Homo_err', 
                'Nonhomo_err', 'Total_err', 
                'Anchor_rt', 'Adj_contig', 'Clusters']
    
    data1 = [
        ['mean/6', '-', '-', '-', '-', '-', '-', '-'],
        ['mean/8', '0.8552', '7.8676', '5.0622', '12.9298', '87.5251', '74.8515', '48'],
        ['mean/10', '0.7736', '9.8756', '10.6296', '20.5052', '82.3276', '63.6886', '48']
    ]
    
    # 表格2数据（进一步缩写）
    headers2 = ['Method', 'Contiguity', 'Homo_err', 
                'Nonhomo_err', 'Total_err', 
                'Anchor_rt', 'Adj_contig']
    
    data2 = [
        ['Random mask', '0.8732', '5.4975', '4.9570', '10.4545', '89.9603', '78.5582'],
        ['w/o Random mask', '0.8559', '9.3436', '8.0841', '17.4277', '89.5845', '76.6754'],
        ['w/o Random mask\n(no mean/6)', '0.8663', '6.0501', '4.8885', '10.9386', '89.6536', '77.6669']
    ]
    
    def create_table_in_ax(ax, headers, data, title, color_index_offset=0):
        """在指定坐标轴上创建表格，使用Paired配色"""
        n_rows = len(data) + 1
        n_cols = len(headers)
        
        # 根据列数调整单元格宽度
        if n_cols == 8:
            cell_width = 0.095
        elif n_cols == 7:
            cell_width = 0.105
        else:
            cell_width = 0.12
        
        cell_height = 0.12  # 增加行高，适应更少的行数
        table_left = 0.05
        table_bottom = 0.1  # 调整底部边距
        
        # 绘制表头 - 使用深色
        header_y = table_bottom + (n_rows - 1) * cell_height
        header_color = colors[7]  # 深紫色作为统一表头颜色
        for col in range(n_cols):
            x = table_left + col * cell_width
            rect = Rectangle((x, header_y), cell_width, cell_height,
                            facecolor=header_color, edgecolor='black', linewidth=1.0)
            ax.add_patch(rect)
            
            # 调整字体大小
            fontsize = 8 if len(headers[col]) > 8 else 9
            ax.text(x + cell_width/2, header_y + cell_height/2, headers[col],
                    ha='center', va='center', color='white', 
                    fontweight='bold', fontsize=fontsize)
        
        # 绘制数据行
        for row_idx, row_data in enumerate(data):
            y_pos = table_bottom + (n_rows - row_idx - 2) * cell_height
            
            # 使用Paired配色交替行颜色
            if row_idx % 2 == 0:
                row_color = colors[0 + color_index_offset]  # 第一组颜色
            else:
                row_color = colors[1 + color_index_offset]  # 第二组颜色
            
            for col in range(n_cols):
                x = table_left + col * cell_width
                rect = Rectangle((x, y_pos), cell_width, cell_height,
                                facecolor=row_color, edgecolor='black', linewidth=0.8,
                                alpha=0.7)
                ax.add_patch(rect)
                
                cell_value = row_data[col]
                
                # 调整字体大小
                if col == 0:
                    fontsize = 8
                    font_weight = 'bold'
                else:
                    fontsize = 8.5
                    font_weight = 'normal'
                
                ax.text(x + cell_width/2, y_pos + cell_height/2, cell_value,
                        ha='center', va='center', color='black', 
                        fontweight=font_weight, fontsize=fontsize)
        
        # 设置坐标轴范围
        table_width = n_cols * cell_width
        table_height = n_rows * cell_height
        ax.set_xlim(table_left - 0.01, table_left + table_width + 0.01)
        ax.set_ylim(table_bottom - 0.01, table_bottom + table_height + 0.1)  # 增加顶部空间
        
        # 添加子标题
        ax.text(0.5, 0.92, title, ha='center', va='top', 
                fontsize=11, fontweight='bold', color='black', transform=ax.transAxes)
    
    # 创建两个表格，使用不同的颜色偏移量
    create_table_in_ax(axes[0], headers1, data1, '(a)', 0)
    create_table_in_ax(axes[1], headers2, data2, '(b)', 2)
    
    plt.tight_layout()
    plt.savefig('combined_tables_summary_paired.png', dpi=300, bbox_inches='tight')
    print("✓ Combined tables saved as 'combined_tables_summary_paired.png' (only (a) and (b))")
    plt.show()

def main():
    """主函数"""
    print("开始创建C88 N50=2M分析表格（Paired配色 + 缩写）...")
    print("="*60)
    
    # 打印Paired配色方案的颜色
    print("Paired配色方案颜色：")
    for i, color in enumerate(colors):
        print(f"  colors[{i}]: {color}")
    print()
    
    # 创建单个表格
    create_table1()
    create_table2()
    create_table3()
    
    # 创建汇总表格
    create_combined_tables()
    
    print("\n" + "="*60)
    print("所有Paired配色缩写表格已生成完成！")
    print("单个表格：")
    print("  - table1_resolution_comparison_paired.png")
    print("  - table2_random_mask_comparison_paired.png")
    print("  - table3_contig_statistics_paired.png")
    print("汇总表格：")
    print("  - combined_tables_summary_paired.png")
    print("="*60)

if __name__ == "__main__":
    main()