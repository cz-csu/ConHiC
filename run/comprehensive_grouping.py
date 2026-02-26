from itertools import product
from collections import defaultdict, OrderedDict
from ortools.sat.python import cp_model
from collections import defaultdict
import argparse
def cluster_with_threshold(result_clusters,fa_dict,cluster_set, unclust_group, unclust_set, size_threshold, similarity_threshold=0.5):
    """
    聚类算法实现
    
    参数:
        cluster_set: 已聚类的集合(set)
        unclust_group: 未聚类的类别名列表(list)
        unclust_set: 对应未聚类的集合(dict: 类别名->元素集合)
        size_threshold: 聚类集合的目标大小
        similarity_threshold: 合并相似度阈值
        
    返回:
        更新后的cluster_set
    """
    
    # 第一阶段：合并未聚类组中相似的类别
    while True:
        max_similarity = -1
        best_pair = None
        
        # 计算所有两两组合的交并比(Jaccard相似度)
        for i in range(len(unclust_group)):
            for j in range(i+1, len(unclust_group)):
                set1 = unclust_set[unclust_group[i]]
                set2 = unclust_set[unclust_group[j]]
                intersection = set1 & set2
                union = set1 | set2
                
                if len(union) == 0:
                    continue
                    
                similarity = len(intersection) / len(union)
                
                if similarity > max_similarity:
                    max_similarity = similarity
                    best_pair = (i, j)
        
        # 如果没有达到阈值的相似对，停止合并
        if max_similarity < similarity_threshold or best_pair is None:
            break
            
        # 合并相似度最高的两个类别
        i, j = best_pair
        new_name = f"merged_{unclust_group[i]}_{unclust_group[j]}"
        new_set = unclust_set[unclust_group[i]] & unclust_set[unclust_group[j]]
        
        # 更新未聚类组和集合
        unclust_group.append(new_name)
        unclust_set[new_name] = new_set
        
        # 删除原来的两个类别
        del unclust_set[unclust_group[i]]
        del unclust_set[unclust_group[j]]
        unclust_group.pop(j)  # 先删除后面的索引
        unclust_group.pop(i)
    print(unclust_set)
    # 第二阶段：按大小顺序将类别与cluster_set合并
    # 按集合大小降序排序
    sorted_groups = sorted(unclust_group, key=lambda x: -len(unclust_set[x]))
    
    for group in sorted_groups:
        if len(cluster_set) >= size_threshold:
            print("Reached size threshold, stopping clustering.")
            break
            
        current_set = unclust_set[group]
        intersection = cluster_set & current_set
        
        # 将交集元素从当前集合中移除
        remaining_elements = current_set - intersection
        if len(remaining_elements) == 0:
            # 如果没有交集，直接跳过
            print(f"remaining_elements for group {group}, skipping.")
            continue
        # 将剩余元素添加到cluster_set
        cluster_set.update(remaining_elements)

        global tot
        result_clusters[tot][0]= remaining_elements
        for frag in remaining_elements:
            result_clusters[tot][1] += fa_dict[frag][1]
        tot+=1
        # 更新未聚类集合(可选)
        #unclust_set[group] = intersection
    print("cluster_set2: ",len(cluster_set))
    return result_clusters
def calculate_total_nContigs(group_to_nContigs):
    return sum(group_to_nContigs.values())
def parse_group_file(file_path):
    """
    解析文件，将 group 映射到 nContigs 和 Contigs。
    
    :param file_path: 文件路径
    :return: 两个字典 (group_to_nContigs, group_to_Contigs)
    """
    group_to_nContigs = {}
    group_to_Contigs = {}

    with open(file_path, "r") as file:
        for line in file:
            if line.startswith("#"):  # 跳过标题行
                continue
            parts = line.strip().split("\t")
            group = parts[0]
            nContigs = int(parts[1])
            contigs = parts[2].split(" ")
            
            group_to_nContigs[group] = nContigs
            group_to_Contigs[group] = contigs

    return group_to_nContigs, group_to_Contigs
def find_dual_intersection_groups(classifications, majority_threshold=0.5):
    """
    找到n个分类集合中交集最大的s个n元组（s是最小分类数），并返回:
    - 最佳n元组列表
    - 未被聚类的组名集合
    - 未被聚类的组内容字典
    
    参数:
        classifications: 字典，键为分类方法名，值为字典{类名: 字符串集合}
        majority_threshold: 多数共识阈值(默认0.5)
        
    返回:
        tuple: (results, unclust_group, unclust_set)
        results: 每个元素是元组(类名1,集合1, 类名2,集合2, ..., 严格交集, 多数共识集合)
        unclust_group: 所有未被选中的类名集合（格式为"方法名_类名"）
        unclust_set: 字典{"方法名_类名": 元素集合}
    """
    n = len(classifications)
    if n == 0:
        return [], set(), {}  # 返回空集合
    
    method_names = list(classifications.keys())
    
    # 确定最小分类数s
    s = min(len(classes) for classes in classifications.values())
    print("ssssssss")
    print(s)
    # 预计算所有类的元素集合
    class_sets = {
        method: {
            cls: set(elements) 
            for cls, elements in classes.items()
        }
        for method, classes in classifications.items()
    }
    
    # 初始化结果和未聚类数据
    results = []
    unclust_group = set()  # 改为集合
    unclust_set = {}
    
    # 已使用的类别，每个分类方法一个集合
    used_classes = {method: set() for method in method_names}
    
    # 找到s个最佳n元组
    for _ in range(s):
        max_intersection_size = -1
        best_combination = None
        best_sets = None
        best_intersection = None
        best_majority_set = None
        
        # 生成候选组合（从第一个分类方法开始）
        base_method = method_names[0]
        
        for base_class in class_sets[base_method]:
            if base_class in used_classes[base_method]:
                continue
                
            current_set = class_sets[base_method][base_class]
            current_combination = [base_class]
            current_sets = [class_sets[base_method][base_class]]
            
            # 为其他分类方法找到最佳匹配类
            for method in method_names[1:]:
                best_class = None
                best_class_set = None
                best_intersection_size = -1
                
                for cls in class_sets[method]:
                    if cls in used_classes[method]:
                        continue
                        
                    intersection = current_set & class_sets[method][cls]
                    if len(intersection) > best_intersection_size:
                        best_intersection_size = len(intersection)
                        best_class = cls
                        best_class_set = class_sets[method][cls]
                
                if best_class is None:
                    break
                    
                current_set = current_set & best_class_set
                current_combination.append(best_class)
                current_sets.append(best_class_set)
            
            # 检查是否找到完整组合
            if len(current_combination) == n:
                # 计算严格交集和多数共识
                intersection = current_set
                intersection_size = len(intersection)
                
                element_counts = defaultdict(int)
                for se in current_sets:
                    for elem in se:
                        element_counts[elem] += 1
                
                majority_set = {
                    elem for elem, count in element_counts.items()
                    if count / n > majority_threshold
                }
                
                if intersection_size > max_intersection_size:
                    max_intersection_size = intersection_size
                    best_combination = current_combination
                    best_sets = current_sets
                    best_intersection = intersection
                    best_majority_set = majority_set
        
        if best_combination is None:
            print("breakbreak")
            break
            
        # 添加到结果
        output_tuple = []
        for cls, se in zip(best_combination, best_sets):
            output_tuple.extend([cls, se])
        output_tuple.append(best_intersection)
        output_tuple.append(best_majority_set)
        results.append(tuple(output_tuple))
        
        # 标记为已使用
        for i, method in enumerate(method_names):
            used_classes[method].add(best_combination[i])
    
    # 收集未聚类的组（使用集合存储）
    for method in method_names:
        for cls in class_sets[method]:
            if cls not in used_classes[method]:
                #full_name = f"{method}_{cls}"
                unclust_group.add(cls)  # 使用add方法添加到集合
                unclust_set[cls] = class_sets[method][cls]
    
    return results, unclust_group, unclust_set
def find_majority_consensus_groups(classifications, majority_threshold=0.5):
    """
    找到n个分类集合中满足多数共识的48个n元组
    
    参数:
        classifications: 字典，键为分类方法名(A1,A2,...)，值为字典{类名: 字符串集合}
        majority_threshold: 多数共识阈值(默认0.5表示超过一半类别包含该元素)
        
    返回:
        列表，每个元素是一个元组，格式为:
        (类名1, 集合1, 类名2, 集合2, ..., 类名n, 集合n, 多数共识集合)
    """
    n = len(classifications)
    if n == 0:
        return []
    
    method_names = list(classifications.keys())
    num_classes = 48  # 每个分类方法有48个类
    results = []
    
    # 已使用的类别，每个分类方法一个集合
    used_classes = {method: set() for method in method_names}
    
    # 预计算所有类的元素集合
    class_sets = {
        method: {
            cls: set(elements) 
            for cls, elements in classes.items()
        }
        for method, classes in classifications.items()
    }
    
    # 找到48个最佳n元组
    for _ in range(num_classes):
        best_combination = None
        best_sets = None
        best_majority_set = None
        best_majority_size = -1
        
        # 生成所有可能的候选组合（避免重复使用已选类别）
        base_method = method_names[0]
        
        for base_class in class_sets[base_method]:
            if base_class in used_classes[base_method]:
                continue
                
            current_combination = [base_class]
            current_sets = [class_sets[base_method][base_class]]
            
            # 为其他分类方法找到最佳匹配类
            for method in method_names[1:]:
                best_class = None
                best_class_set = None
                best_majority = -1
                
                for cls in class_sets[method]:
                    if cls in used_classes[method]:
                        continue
                        
                    # 临时组合当前选择的类
                    temp_combination = current_combination + [cls]
                    temp_sets = current_sets + [class_sets[method][cls]]
                    
                    # 计算多数共识集合
                    element_counts = defaultdict(int)
                    for s in temp_sets:
                        for elem in s:
                            element_counts[elem] += 1
                    
                    majority_set = {
                        elem for elem, count in element_counts.items()
                        if count / len(temp_sets) > majority_threshold
                    }
                    
                    majority_size = len(majority_set)
                    
                    if majority_size > best_majority:
                        best_majority = majority_size
                        best_class = cls
                        best_class_set = class_sets[method][cls]
                
                if best_class is None:
                    break
                    
                current_combination.append(best_class)
                current_sets.append(best_class_set)
            
            # 如果成功找到完整组合
            if len(current_combination) == n:
                # 计算最终多数共识集合
                element_counts = defaultdict(int)
                for s in current_sets:
                    for elem in s:
                        element_counts[elem] += 1
                
                majority_set = {
                    elem for elem, count in element_counts.items()
                    if count / n > majority_threshold
                }
                
                majority_size = len(majority_set)
                
                if majority_size > best_majority_size:
                    best_majority_size = majority_size
                    best_combination = current_combination
                    best_sets = current_sets
                    best_majority_set = majority_set
        
        if best_combination is None:
            break
            
        # 构建输出元组：交替类名和集合，最后加多数共识集合
        output_tuple = []
        for cls, s in zip(best_combination, best_sets):
            output_tuple.append(cls)
            output_tuple.append(s)
        output_tuple.append(best_majority_set)
        
        results.append(tuple(output_tuple))
        
        # 标记这些类为已使用
        for i, method in enumerate(method_names):
            used_classes[method].add(best_combination[i])
    
    return results   
def find_optimal_triplets(contigs, n=48):
    # 生成所有可能的三元组 (A_i, B_j, C_k)
    triplets = []
    n_methods = len(contigs)
    for k in range(n_methods):
        triplet = []
        for index, (ai, a_set) in enumerate(contigs[k].items()):
            if index == 0:  # 跳过第一个 item
                common= set(a_set)
            else:
                common = common & set(a_set)
            triplet.append(ai)
            triplet.append(a_set)
        #triplets.append((len(common), ai, bj, ck ,a_set ,b_set ,c_set ,common))
        triplet.append(common)
        triplet.append(len(common))
        triplets.append(triplet)
    # 按交集大小降序排序
    triplets.sort(reverse=True, key=lambda x: x[-1])
    
    selected = [set() for _ in range(n_methods)]
    result = []
    
    for triplet in triplets:
        flag=0
        for i in range(0, len(triplets)-2, 2):
            if triplet[i] in selected[i//2]:
                flag=1
                break
        if flag==0:
            for i in range(0, len(triplets)-2, 2):
                selected[i//2].add(triplet[i]) 
            result.append(triplet)
            if len(result) == n:
                break
    return result
def solve_3d_assignment_ortools(A, B, C, n=48):
    model = cp_model.CpModel()

    # 变量：x[i][j][k] = 1 表示选择 (A_i, B_j, C_k)
    x = {}
    for i in A:
        for j in B:
            for k in C:
                x[(i, j, k)] = model.NewBoolVar(f"x_{i}_{j}_{k}")

    # 目标函数：最大化交集大小之和
    model.Maximize(
        sum(len(A[i] & B[j] & C[k]) * x[(i, j, k)] for i in A for j in B for k in C)
    )

    # 约束1：每个 A_i 只能选一次
    for i in A:
        model.Add(sum(x[(i, j, k)] for j in B for k in C) <= 1)

    # 约束2：每个 B_j 只能选一次
    for j in B:
        model.Add(sum(x[(i, j, k)] for i in A for k in C) <= 1)

    # 约束3：每个 C_k 只能选一次
    for k in C:
        model.Add(sum(x[(i, j, k)] for i in A for j in B) <= 1)

    # 约束4：总共选 n 个三元组
    model.Add(sum(x[(i, j, k)] for i in A for j in B for k in C) == n)

    # 求解
    solver = cp_model.CpSolver()
    status = solver.Solve(model)

    # 提取结果
    solution = []
    if status == cp_model.OPTIMAL:
        for i in A:
            for j in B:
                for k in C:
                    if solver.Value(x[(i, j, k)]) == 1:
                        common = len(A[i] & B[j] & C[k])
                        solution.append((i, j, k, common))
    return solution
def write_group_data_to_file(group_to_nContigs, group_to_Contigs, output_file):
    """
    将两个字典的数据写入到文件中，格式为 #Group nContigs Contigs。
    
    :param group_to_nContigs: 字典，key 为 group，value 为 nContigs
    :param group_to_Contigs: 字典，key 为 group，value 为 Contigs 列表
    :param output_file: 输出文件路径
    """
    with open(output_file, "w") as file:
        # 写入标题行
        file.write("#Group\tnContigs\tContigs\n")
        
        # 遍历 group_to_nContigs 的 key
        for group in group_to_nContigs:
            nContigs = group_to_nContigs[group]
            contigs = " ".join(group_to_Contigs[group])  # 将 Contigs 列表拼接成字符串
            file.write(f"{group}\t{nContigs}\t{contigs}\n")
def parse_RE_sites(sites):

    output_sites = list()

    for site in sites:
        if 'N' in site:
            output_sites.append(site.replace('N', 'A', 1))
            output_sites.append(site.replace('N', 'T', 1))
            output_sites.append(site.replace('N', 'C', 1))
            output_sites.append(site.replace('N', 'G', 1))
        else:
            output_sites.append(site)

    if 'N' not in ''.join(output_sites):
        return output_sites
    else:
        return parse_RE_sites(output_sites)
def count_RE_sites(seq, RE):

    sites = [site.strip().upper() for site in RE.split(',') if site.strip()]
    parsed_sites = parse_RE_sites(sites)

    RE_sites = 0
    for site in parsed_sites:
        RE_sites += seq.count(site)

    return RE_sites
def parse_fasta(fasta, RE='GATC', keep_letter_case=False):

    fa_dict = dict()
    with open(fasta) as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith('>'):
                ctg = line.split()[0][1:]
                fa_dict[ctg] = list()
            else:
                if keep_letter_case:
                    fa_dict[ctg].append(line.strip())
                else:
                    fa_dict[ctg].append(line.strip().upper())

    for ctg, seq_list in fa_dict.items():
        # joining list is faster than concatenating strings
        seq = ''.join(seq_list)
        # add pseudo-count of 1 to prevent division by zero (as what ALLHiC does)
        RE_sites = count_RE_sites(seq, RE) + 1
        fa_dict[ctg] = [seq, len(seq), RE_sites]

    return fa_dict
def parse_arguments():

    parser = argparse.ArgumentParser(prog='haphic cluster')

    # Parameters for parsing input files and pipeline control
    input_group = parser.add_argument_group('>>> Parameters for parsing input files and pipeline control')
    input_group.add_argument(
            'data_paths', type=list, default=[], help='data_paths')
    
    args = parser.parse_args()

    return args
def run(data_paths):
    # 调用示例
    #data_paths =["/home/chenzh/HapHiC/data/CS/RM/02.reassign/99218/final_groups/final_clusters.txt",
    #       "/home/chenzh/HapHiC/data/CS/RM/02.reassign/148827/final_groups/final_clusters.txt",
    #        "/home/chenzh/HapHiC/data/CS/RM/02.reassign/198436/final_groups/final_clusters.txt",]
    n_methods=len(data_paths)
    all_contig= set()
    nContigs = [0] * n_methods
    Contigs = [0] * n_methods
    total_nContigs = [0] * n_methods
    size_threshold = 0
    for i in range(n_methods):
        nContigs[i], Contigs[i] = parse_group_file(data_paths[i])
        for ctg in Contigs[i].values():
            all_contig |= set(ctg)
        total_nContigs[i] = calculate_total_nContigs(nContigs[i])
        #size_threshold=max(size_threshold, total_nContigs[i])
    print("all_contig: ", len(all_contig))
    size_threshold=len(all_contig) 
    print("size_threshold: ",size_threshold)
    #total_A_nContigs = calculate_total_nContigs(A_nContigs)
    #total_B_nContigs = calculate_total_nContigs(B_nContigs)
    #total_C_nContigs = calculate_total_nContigs(C_nContigs)
    #print(total_A_nContigs,total_B_nContigs,total_C_nContigs)
    #print("Group to nContigs:", group_to_nContigs)
    #print("Group to Contigs:", group_to_Contigs)
    # 调用函数
    Contigs_dict = {index: value for index, value in enumerate(Contigs)}
    optimal_triplets, unclust_group, unclust_set = find_dual_intersection_groups(Contigs_dict)
    #optimal_triplets = find_majority_consensus_groups(Contigs_dict)
    """
    A_Contigs = {k: set(v) for k, v in A_Contigs.items()}
    B_Contigs = {k: set(v) for k, v in B_Contigs.items()}
    C_Contigs = {k: set(v) for k, v in C_Contigs.items()}
    solution = solve_3d_assignment_ortools(A_Contigs, B_Contigs, C_Contigs)
    for triplet in solution:
        print(triplet)
    """
    sum_num=0
    group_to_nContigs={}
    group_to_Contigs={}
    output_file = "output.txt"
    fasta_file = "/home/chenzh/HapHiC/data/ap/shuffled_ap_split_genome.fa"
    outdir = "."
    print("000000000000000000000")
    print(len(optimal_triplets))
    for n,triplet in enumerate(optimal_triplets):
        print(len(triplet))
        odd_position_values=[]
        for index, value in enumerate(triplet):
            if index==len(triplet)-1:
                odd_position_values.append(len(value))
                break                
            if index==len(triplet)-2:
                odd_position_values.append(len(value))
                sum_num+=len(value)
                continue
            if index % 2 == 0:
                odd_position_values.append(value)
            else:
                odd_position_values.append(len(value))
        #odd_position_values = [value for index, value in enumerate(triplet) if index % 2 == 0]
        print(odd_position_values)
    print("sum_num: ", sum_num) 
    #"""
    fa_dict = parse_fasta(fasta_file,RE="AAGCTT")
    print(len(fa_dict))
    cluster_set = set()
    #unclust_group = set()
    #unclust_set = {}
    #print(fa_dict)
    result_clusters = defaultdict(lambda: [[], 0])
    global tot
    tot=0
    adapt_set=-1 #-1是多数共识 -2是严格交集
    for n,triplet in enumerate(optimal_triplets):
        if len(triplet[adapt_set])!=0:
            result_clusters[tot][0]= triplet[adapt_set]
            for frag in triplet[adapt_set]:
                result_clusters[tot][1] += fa_dict[frag][1]
            tot+=1
            cluster_set|=set(triplet[adapt_set])
        else:
            for index, value in enumerate(triplet):
                if index==len(triplet)-2:
                    break
                if index % 2 == 0:
                    unclust_group.add(value)
                    unclust_set[value]=set(triplet[index+1])
    print("cluster_set: ",len(cluster_set))
    print(tot)
    print("unclust_group_len: ",len(unclust_group))
    unclust_group=list(unclust_group)
    # 执行聚类
    print("开始聚类... ",len(result_clusters))
    print(len(unclust_group))
    result_clusters = cluster_with_threshold(
        result_clusters,
        fa_dict,
        cluster_set.copy(),  # 使用副本以避免修改原集合
        unclust_group.copy(),
        {k: v.copy() for k, v in unclust_set.items()},  # 深拷贝
        size_threshold=size_threshold,
        similarity_threshold=0.5
    )
    result_clusters = sorted(tuple(result_clusters.values()), key=lambda x: x[1], reverse=True)
    print("最终聚类结果数: ", len(result_clusters))
    #print("最终聚类结果:", result)
        #print(triplet)
    # output a cluster file for result overview and reassignment
    #"""
    tot_num=0
    with open('{0}/mcl_{0}.clusters.txt'.format(outdir), 'w') as fout:
        fout.write('#Group\tnContigs\tContigs\n')
        for n, (ctgs, group_len) in enumerate(result_clusters, 1):
            # sort contigs by length in each cluster
            result_clusters[n-1][0]=list(result_clusters[n-1][0])
            result_clusters[n-1][0].sort(key=lambda x: fa_dict[x][1], reverse=True)
            fout.write('group{}_{}bp\t{}\t{}\n'.format(n, group_len, len(ctgs), ' '.join(ctgs)))
            tot_num+=len(ctgs)
            # generate group*.txt files for ALLHiC optimize
    print("Final num: ",tot_num)
            # generate group*.txt files for ALLHiC optimize
    """
    for n, (ctgs, group_len) in enumerate(result_clusters, 1):
        with open('{}/group{}_{}bp.txt'.format(outdir, n, group_len), 'w') as fout:
            fout.write('#Contig\tRECounts\tLength\n')
            for ctg in ctgs:
                length, RE_sites = fa_dict[ctg][1:3]
                fout.write('{}\t{}\t{}\n'.format(ctg, RE_sites, length))
    #write_group_data_to_file(group_to_nContigs, group_to_Contigs, output_file)
    #"""
    return '{0}/mcl_{0}.clusters.txt'.format(outdir)
def main():

    # get arguments
    #args = parse_arguments()
    #run(args, 'HapHiC_cluster.log')
    run([])  # 传入空列表作为示例

if __name__ == "__main__":

    main()
