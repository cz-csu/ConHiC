import os
from collections import Counter, defaultdict
from tqdm import tqdm
import numpy as np

def process_group_files(files, output_dir, name, base_file_index=0):
    """处理一组文件，基于base文件对齐，选择平均边支持度最高的文件"""
    
    if not files:
        print(f"No files for {name}")
        return
    
    print(f"\nProcessing {name}")
    print(f"Number of files: {len(files)}")
    
    # 存储所有对齐后的文件内容
    aligned_contents = []
    file_paths = []
    
    # 1. 读取base文件
    if base_file_index >= len(files):
        print(f"  Warning: base_file_index {base_file_index} out of range, using 0")
        base_file_index = 0
    
    base_file = files[base_file_index]
    try:
        with open(base_file, 'r') as f:
            # 读取文件的最后一行（最后一个非空行）
            last_line = None
            for line in f:
                line = line.strip()
                if line:  # 只保存非空行
                    last_line = line
            
            if not last_line:
                print(f"  Error: Base file {base_file} has no valid content")
                return
            
            base_elements = last_line.split()
        
        if len(base_elements) < 2:
            print(f"  Error: Base file has too few elements")
            return
        
        # 存储base文件（不对其进行方向反转）
        aligned_contents.append({
            'idx': base_file_index,
            'original_idx': base_file_index,
            'path': base_elements,
            'file_path': base_file,
            'elements': set([elem[:-1] if elem[-1] in '+-' else elem 
                           for elem in base_elements]),
            'is_reversed': False,
            'similarity': 1.0  # base文件与自身的相似度为1
        })
        file_paths.append(base_file)
        
        print(f"Base file: index {base_file_index}, length: {len(base_elements)}")
        
    except Exception as e:
        print(f"  Error reading base file {base_file}: {str(e)}")
        return
    
    # 提取base文件的元素（不带方向）
    base_elements_no_dir = [elem[:-1] if elem[-1] in '+-' else elem 
                           for elem in base_elements]
    
    # 2. 处理其他文件，与base对齐
    for file_idx, file in enumerate(files):
        if file_idx == base_file_index:
            continue  # base文件已处理
        try:
            with open(file, 'r') as f:
                # 读取文件的最后一行（最后一个非空行）
                last_line = None
                for line in f:
                    line = line.strip()
                    if line:  # 只保存非空行
                        last_line = line
                
                if not last_line:
                    print(f"  Warning: No valid content in {file}")
                    continue
                
                elements_with_direction = last_line.split()
            
            if len(elements_with_direction) < 2:
                print(f"  Warning: Too few elements in {file}")
                continue
            
            # 提取元素（不带方向）
            elements_no_dir = [elem[:-1] if elem[-1] in '+-' else elem 
                             for elem in elements_with_direction]
            
            # 计算与base文件的相似度（正反两个方向）
            # 获取公共元素
            common_elements = set(base_elements_no_dir) & set(elements_no_dir)
            
            if len(common_elements) < 2:
                print(f"  File {file_idx}: Not enough common elements with base ({len(common_elements)})")
                # 如果公共元素太少，直接使用原方向
                aligned_path = elements_with_direction
                is_reversed = False
                similarity = len(common_elements) / max(len(base_elements_no_dir), len(elements_no_dir))
            else:
                # 创建base文件中公共元素的位置映射
                base_order = {elem: idx for idx, elem in enumerate(base_elements_no_dir) 
                             if elem in common_elements}
                
                # 正向序列中公共元素的位置映射
                forward_order = {elem: idx for idx, elem in enumerate(elements_no_dir) 
                                if elem in common_elements}
                
                # 计算正向相似度（顺序一致的对数比例）
                forward_sim = sum(1 for elem1 in common_elements for elem2 in common_elements 
                                if elem1 != elem2 and 
                                (base_order[elem1] < base_order[elem2]) == (forward_order[elem1] < forward_order[elem2]))
                
                # 总可能的比较对数
                total_pairs = len(common_elements) * (len(common_elements) - 1)
                forward_similarity = forward_sim / total_pairs if total_pairs > 0 else 0
                
                # 计算反向相似度
                reversed_elements = elements_with_direction[::-1]
                reversed_elements_no_dir = elements_no_dir[::-1]
                reversed_order = {elem: idx for idx, elem in enumerate(reversed_elements_no_dir) 
                                if elem in common_elements}
                
                reverse_sim = sum(1 for elem1 in common_elements for elem2 in common_elements 
                                if elem1 != elem2 and 
                                (base_order[elem1] < base_order[elem2]) == (reversed_order[elem1] < reversed_order[elem2]))
                reverse_similarity = reverse_sim / total_pairs if total_pairs > 0 else 0
                
                # 根据相似度决定是否反转
                if reverse_similarity > forward_similarity:
                    # 反转序列并反转方向
                    reversed_with_inverted_dir = []
                    for elem in reversed_elements:
                        if elem[-1] == '+':
                            reversed_with_inverted_dir.append(elem[:-1] + '-')
                        elif elem[-1] == '-':
                            reversed_with_inverted_dir.append(elem[:-1] + '+')
                        else:
                            reversed_with_inverted_dir.append(elem + '-')  # 默认
                    aligned_path = reversed_with_inverted_dir
                    is_reversed = True
                    similarity = reverse_similarity
                else:
                    aligned_path = elements_with_direction
                    is_reversed = False
                    similarity = forward_similarity
            
            # 存储对齐后的内容
            aligned_contents.append({
                'idx': len(aligned_contents),  # 新索引
                'original_idx': file_idx,
                'path': aligned_path,
                'file_path': file,
                'elements': set([elem[:-1] if elem[-1] in '+-' else elem 
                               for elem in aligned_path]),
                'is_reversed': is_reversed,
                'similarity': similarity
            })
            file_paths.append(file)
            
            print(f"  File {file_idx}: aligned with base, similarity={similarity:.3f}, "
                  f"reversed={is_reversed}, length={len(aligned_path)}")
            
        except Exception as e:
            print(f"  Error processing file {file}: {str(e)}")
            continue
    
    if len(aligned_contents) < 1:
        print(f"  Error: No valid aligned files")
        return
    
    print(f"\nSuccessfully aligned {len(aligned_contents)} files")
    
    # 3. 计算每个文件的边支持度
    print("\nCalculating edge support...")
    
    # 首先收集所有文件的所有边
    all_edges_lists = []
    for content in aligned_contents:
        path = content['path']
        edges = []
        for i in range(len(path)-1):
            edges.append((path[i], path[i+1]))
        all_edges_lists.append(edges)
    
    # 统计所有边在所有对齐文件中的出现次数
    edge_counter = Counter()
    for edges in all_edges_lists:
        for edge in edges:
            edge_counter[edge] += 1
    
    print(f"Total unique edges across aligned files: {len(edge_counter)}")
    
    # 4. 为每个对齐后的文件计算平均边支持度
    file_scores = []
    for i, content in enumerate(aligned_contents):
        edges = all_edges_lists[i]
        
        if len(edges) == 0:
            avg_score = 0
        else:
            edge_scores = [edge_counter[edge] for edge in edges]
            avg_score = np.mean(edge_scores)
        
        # 计算覆盖率（与base文件的公共元素比例）
        coverage = len(content['elements'] & aligned_contents[0]['elements']) / len(aligned_contents[0]['elements'])
        
        file_scores.append({
            'idx': content['idx'],
            'original_idx': content['original_idx'],
            'score': avg_score,
            'coverage': coverage,
            'similarity': content['similarity'],
            'edge_count': len(edges),
            'path_length': len(content['path']),
            'unique_elements': len(content['elements']),
            'is_reversed': content['is_reversed'],
            'is_base': (i == 0)
        })
    
    # 5. 按得分排序，选择得分最高的文件
    # 排序规则：1) 平均边支持度 2) 覆盖率 3) 与base的相似度
    file_scores.sort(key=lambda x: (-x['score'], -x['coverage'], -x['similarity']))
    
    # 打印所有文件的得分
    print("\nAligned file scores (sorted by average edge support):")
    for i, score_info in enumerate(file_scores):
        base_marker = " [BASE]" if score_info['is_base'] else ""
        reverse_marker = " [REV]" if score_info['is_reversed'] else ""
        print(f"  {i+1}. File {score_info['original_idx']}{base_marker}{reverse_marker}: "
              f"avg_edge={score_info['score']:.3f}, "
              f"coverage={score_info['coverage']:.3f}, "
              f"similarity={score_info['similarity']:.3f}, "
              f"elements={score_info['unique_elements']}")
    
    # 选择最佳文件
    best_score_info = file_scores[0]
    best_file_idx = best_score_info['idx']
    best_file = next(content for content in aligned_contents if content['idx'] == best_file_idx)
    best_path = best_file['path']
    
    print(f"\nSelected file: original index {best_score_info['original_idx']}")
    if best_score_info['is_base']:
        print("  (Base file selected)")
    if best_score_info['is_reversed']:
        print("  (File was reversed during alignment)")
    print(f"  Average edge support: {best_score_info['score']:.3f}")
    print(f"  Coverage with base: {best_score_info['coverage']:.3f}")
    print(f"  Similarity with base: {best_score_info['similarity']:.3f}")
    print(f"  Path length: {len(best_path)}")
    print(f"  Unique elements: {best_score_info['unique_elements']}")
    
    # 6. 计算最终路径的详细边支持度统计
    print("\nFinal path edge support distribution:")
    final_edges = []
    final_edge_supports = []
    
    for i in range(len(best_path)-1):
        edge = (best_path[i], best_path[i+1])
        final_edges.append(edge)
        support = edge_counter.get(edge, 0)
        final_edge_supports.append(support)
    
    if final_edge_supports:
        # 计算支持度分布
        support_counts = Counter(final_edge_supports)
        total_edges = len(final_edges)
        
        print(f"  Total edges: {total_edges}")
        for support_level in sorted(support_counts.keys()):
            count = support_counts[support_level]
            percentage = count / total_edges * 100
            print(f"    {support_level}-support edges: {count} ({percentage:.1f}%)")
        
        # 计算统计量
        avg_support = np.mean(final_edge_supports)
        min_support = np.min(final_edge_supports)
        max_support = np.max(final_edge_supports)
        std_support = np.std(final_edge_supports)
        
        print(f"\n  Statistics:")
        print(f"    Average support: {avg_support:.2f}")
        print(f"    Minimum support: {min_support}")
        print(f"    Maximum support: {max_support}")
        print(f"    Standard deviation: {std_support:.2f}")
        
        # 完全支持的边（在所有对齐文件中都出现）
        fully_supported = sum(1 for support in final_edge_supports if support == len(aligned_contents))
        if len(aligned_contents) > 0:
            fully_supported_ratio = fully_supported / total_edges if total_edges > 0 else 0
            print(f"    Fully supported edges: {fully_supported}/{total_edges} ({fully_supported_ratio:.1%})")
    
    # 7. 保存结果
    output_file = os.path.join(output_dir, name)
    with open(output_file, 'w') as f:
        f.write(">INIT\n")
        f.write(" ".join(best_path) + "\n")
    
    print(f"\nResults saved to: {output_file}")
    
    # 8. 额外信息：如果最佳文件不是base，显示与base的差异
    if not best_score_info['is_base'] and aligned_contents[0]['path']:
        base_path = aligned_contents[0]['path']
        if len(best_path) == len(base_path):
            # 计算方向差异
            dir_differences = 0
            for i in range(len(best_path)):
                best_elem = best_path[i]
                base_elem = base_path[i] if i < len(base_path) else ""
                
                if base_elem and best_elem:
                    # 提取元素名和方向
                    best_name = best_elem[:-1] if best_elem[-1] in '+-' else best_elem
                    best_dir = best_elem[-1] if best_elem[-1] in '+-' else '+'
                    base_name = base_elem[:-1] if base_elem[-1] in '+-' else base_elem
                    base_dir = base_elem[-1] if base_elem[-1] in '+-' else '+'
                    
                    if best_name == base_name and best_dir != base_dir:
                        dir_differences += 1
            
            if dir_differences > 0:
                print(f"  Note: Selected path has {dir_differences} elements with different direction from base")

def main():
    k = 4
    input_dirs = [f"/home/chenzh/HapHiC/data/ap/RM/oursort/fastsort_allhic/1.0{i}/final_tours" for i in range(k)]
    output_dir = "/home/chenzh/HapHiC/data/ap/RM/oursort/aggregated_tours"
    
    os.makedirs(output_dir, exist_ok=True)
    
    # 检查输入目录
    print("Checking input directories:")
    valid_dirs = []
    for i, d in enumerate(input_dirs):
        if os.path.exists(d):
            print(f"  Directory {i}: {d} - EXISTS")
            valid_dirs.append(d)
        else:
            print(f"  Directory {i}: {d} - NOT FOUND")
    
    if not valid_dirs:
        print("No valid input directories found!")
        return
    
    # 收集所有文件
    group_files_dict = defaultdict(list)
    for d in valid_dirs:
        if os.path.exists(d):
            for f in sorted(os.listdir(d)):
                if f.endswith('.tour'):
                    group_files_dict[f].append(os.path.join(d, f))
    
    print(f"\nFound {len(group_files_dict)} unique tour files")
    
    # 统计每个group的文件数
    file_counts = Counter()
    for name, files in group_files_dict.items():
        file_counts[len(files)] += 1
    
    print("\nFile count distribution per group:")
    for count in sorted(file_counts.keys()):
        print(f"  {count} files: {file_counts[count]} groups")
    
    # 处理每组文件
    total_groups = len(group_files_dict)
    processed_groups = 0
    
    for name, files in tqdm(group_files_dict.items(), desc="Processing tour files"):
        print(f"\n{'='*60}")
        print(f"Group {processed_groups + 1}/{total_groups}: {name}")
        print(f"{'='*60}")
        
        try:
            process_group_files(files, output_dir, name, base_file_index=0)
            processed_groups += 1
        except Exception as e:
            print(f"\nFailed on {name}: {str(e)}")
            import traceback
            traceback.print_exc()
            continue
    
    print(f"\n{'='*60}")
    print("Selection complete!")
    print(f"{'='*60}")
    print(f"Total groups processed: {processed_groups}/{total_groups}")
    print(f"Output directory: {output_dir}")
    
    # 检查输出文件
    output_files = sorted([f for f in os.listdir(output_dir) if f.endswith('.tour')])
    print(f"\nGenerated {len(output_files)} output files in {output_dir}")

if __name__ == "__main__":
    main()