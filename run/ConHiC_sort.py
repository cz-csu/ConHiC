import os
import sys
import subprocess
import logging
import time
import gc
import shutil

from ConHiC_build import parse_tours

import pickle
import argparse
from collections import defaultdict
from multiprocessing import Pool
from itertools import combinations, product

from numpy import float32, zeros, hstack
from networkx import Graph, connected_components, shortest_path
from networkx import tree as nxtree
from scipy.sparse import coo_matrix
import numpy as np
import warnings
warnings.filterwarnings('ignore', category=FutureWarning)

from _version import __version__, __update_time__

logging.basicConfig(
        format='%(asctime)s <%(filename)s> [%(funcName)s] %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
        )

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def load_sequence_data(fasta_path):
    """
    Extract sequence length information from FASTA file.
    
    Parameters:
    -----------
    fasta_path : str
        Path to the FASTA format genome file
        
    Returns:
    --------
    dict
        Dictionary mapping sequence identifiers to their lengths
    """
    logger.info('Reading sequence file and calculating lengths...')

    seq_lengths = dict()

    with open(fasta_path) as f:
        current_seq = None
        for line in f:
            if not line.strip():
                continue
            if line.startswith('>'):
                current_seq = line.split()[0][1:]
                seq_lengths[current_seq] = 0
            else:
                seq_lengths[current_seq] += len(line.strip())
    return seq_lengths


def sparse_matrix_converter(edge_weights, matrix_dimension, include_self_connections=False):
    """
    Convert edge weight dictionary to dense matrix format.
    
    Parameters:
    -----------
    edge_weights : dict
        Dictionary with tuple keys (i,j) representing edges and float values as weights
    matrix_dimension : int
        Size of the square output matrix
    include_self_connections : bool
        Whether to add self-loops (diagonal elements)
        
    Returns:
    --------
    numpy.ndarray
        Dense matrix representation of the graph
    """
    rows, cols, values = list(), list(), list()

    for (node_i, node_j), weight in edge_weights.items():
        rows.append(node_i)
        cols.append(node_j)
        # Maintain symmetry
        rows.append(node_j)
        cols.append(node_i)
        values.append(weight)
        values.append(weight)

    # Add self-connections if requested
    if include_self_connections:
        for node in range(matrix_dimension):
            rows.append(node)
            cols.append(node)
            values.append(1)

    # Build sparse matrix then convert to dense
    return coo_matrix((values, (rows, cols)), 
                      shape=(matrix_dimension, matrix_dimension), 
                      dtype=float32).toarray()


def process_group_input(group_file, contact_map_dir):
    """
    Parse group file and locate corresponding contact map file.
    
    Parameters:
    -----------
    group_file : str
        Path to group file containing sequence information
    contact_map_dir : str
        Directory containing contact map files (.clm)
        
    Returns:
    --------
    tuple
        (sequence_info_list, contact_map_path, group_prefix)
    """
    group_prefix = os.path.splitext(os.path.basename(group_file))[0]
    contact_map_path = '{}/{}.clm'.format(contact_map_dir, group_prefix)

    if not os.path.exists(contact_map_path):
        raise IOError('Contact map file verification failed: Cannot locate corresponding .clm file '
                'in {} for group file {}'.format(contact_map_dir, group_file))

    sequence_info = list()

    with open(group_file) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.split()
            sequence_info.append((fields[0], int(fields[2])))

    # Sort by sequence length (descending)
    sequence_info.sort(key=lambda x: x[1], reverse=True)

    return sequence_info, contact_map_path, group_prefix


def extract_subgraph_contacts(sequence_list, full_contact_dictionary):
    """
    Extract contact information for a subset of sequences.
    
    Parameters:
    -----------
    sequence_list : list
        List of sequence identifiers
    full_contact_dictionary : dict
        Complete contact dictionary from pickle file
        
    Returns:
    --------
    tuple
        (subgraph_contacts, endpoint_index_mapping)
    """
    endpoint_tags = ['_H', '_T']
    tag_combinations = list(product(endpoint_tags, repeat=2))

    current_index = 0
    subgraph_contacts = defaultdict(int)
    endpoint_index_map = dict()

    for seq_a, seq_b in combinations(sequence_list, 2):
        seq_a, seq_b = sorted([seq_a, seq_b])
        for tag_a, tag_b in tag_combinations:
            endpoint_a, endpoint_b = seq_a + tag_a, seq_b + tag_b
            if (endpoint_a, endpoint_b) in full_contact_dictionary:
                contact_count = full_contact_dictionary[(endpoint_a, endpoint_b)]
            else:
                contact_count = 0
            
            for endpoint in (endpoint_a, endpoint_b):
                if endpoint not in endpoint_index_map:
                    endpoint_index_map[endpoint] = current_index
                    current_index += 1
                    
            if contact_count:
                subgraph_contacts[(endpoint_index_map[endpoint_a], 
                                   endpoint_index_map[endpoint_b])] = contact_count

    return subgraph_contacts, endpoint_index_map


def calculate_endpoint_length(endpoint, sequence_lengths):
    """
    Calculate length of an endpoint (half-sequence for contigs, sum for scaffolds).
    
    Parameters:
    -----------
    endpoint : str or tuple
        Endpoint identifier or scaffold tuple
    sequence_lengths : dict
        Dictionary of sequence lengths
        
    Returns:
    --------
    float
        Length of the endpoint
    """
    if isinstance(endpoint, str):
        return sequence_lengths[endpoint.rsplit('_', 1)[0]] / 2
    else:
        assert isinstance(endpoint, tuple)
        return sum([calculate_endpoint_length(ep, sequence_lengths) for ep in endpoint])


def compute_density_matrix(contact_matrix, matrix_dimension, index_to_endpoint_dict, 
                          sequence_lengths, flanking_dict, density_method):
    """
    Calculate normalized contact density matrix.
    
    Parameters:
    -----------
    contact_matrix : numpy.ndarray
        Raw contact count matrix
    matrix_dimension : int
        Size of the matrix
    index_to_endpoint_dict : dict
        Mapping from indices to endpoint identifiers
    sequence_lengths : dict
        Sequence length information
    flanking_dict : dict
        Flanking region information for scaffolds
    density_method : str
        Method for density calculation ('sum', 'multiplication', or 'geometric_mean')
        
    Returns:
    --------
    numpy.ndarray
        Normalized density matrix
    """
    def get_endpoint_length(endpoint_id):
        if endpoint_id in flanking_dict:
            return flanking_dict[endpoint_id][1]
        else:
            return calculate_endpoint_length(endpoint_id, sequence_lengths)

    length_product_dict = dict()
    indices = index_to_endpoint_dict.keys()

    for idx_i, idx_j in combinations(indices, 2):
        endpoint_i = index_to_endpoint_dict[idx_i]
        endpoint_j = index_to_endpoint_dict[idx_j]
        len_i = get_endpoint_length(endpoint_i)
        len_j = get_endpoint_length(endpoint_j)

        if density_method == 'sum':
            combined_length = len_i + len_j
        elif density_method == 'multiplication':
            combined_length = len_i * len_j
        elif density_method == 'geometric_mean':
            combined_length = (len_i * len_j) ** 0.5

        length_product_dict[(idx_i, idx_j)] = combined_length

    # Add self-loops to avoid division by zero
    length_matrix = sparse_matrix_converter(length_product_dict, matrix_dimension, 
                                            include_self_connections=True)

    return contact_matrix / length_matrix


def calculate_edge_confidence(matrix_dimension, active_edges, density_matrix):
    """
    Calculate confidence scores for graph edges.
    
    Parameters:
    -----------
    matrix_dimension : int
        Size of the matrix
    active_edges : dict
        Dictionary of edges with non-zero weights
    density_matrix : numpy.ndarray
        Normalized density matrix
        
    Returns:
    --------
    tuple
        (confidence_matrix, max_non_sister_density)
    """
    confidence_matrix = zeros((matrix_dimension, matrix_dimension))

    for idx_i, idx_j in active_edges:
        density = density_matrix[idx_i, idx_j]

        # Collect all densities incident on either node (excluding current edge)
        all_densities = hstack((density_matrix[idx_i, :idx_j], 
                               density_matrix[idx_i, idx_j+1:], 
                               density_matrix[:, idx_j]))
        
        # Find maximum and second maximum
        max_pos = all_densities.argmax()
        all_densities = hstack((all_densities[:max_pos], all_densities[max_pos+1:]))
        second_max = all_densities.max()
        
        # Calculate confidence score
        if density == 0:
            confidence = 0
        elif second_max == 0:
            confidence = 2
        else:
            confidence = density / second_max
            
        confidence_matrix[idx_i, idx_j] = confidence
        confidence_matrix[idx_j, idx_i] = confidence

    max_non_sister = confidence_matrix.max()

    return confidence_matrix, max_non_sister


def apply_confidence_filter(confidence_matrix, active_edges, threshold):
    """
    Filter edges based on confidence threshold.
    
    Parameters:
    -----------
    confidence_matrix : numpy.ndarray
        Matrix of confidence scores
    active_edges : dict
        Dictionary of active edges
    threshold : float
        Confidence threshold for edge retention
    """
    for idx_i, idx_j in active_edges:
        if confidence_matrix[idx_i, idx_j] <= threshold:
            confidence_matrix[idx_i, idx_j] = 0
            confidence_matrix[idx_j, idx_i] = 0


def extract_endpoint_components(endpoint):
    """
    Extract individual components from an endpoint (contig or scaffold).
    
    Parameters:
    -----------
    endpoint : str or tuple
        Endpoint identifier
        
    Returns:
    --------
    tuple
        (component_tuple, leftmost_component, rightmost_component)
    """
    if isinstance(endpoint, str):
        return (endpoint,), endpoint, endpoint
    else:
        assert isinstance(endpoint, tuple)
        return endpoint, endpoint[0], endpoint[-1]


def partition_scaffold_path(path, sequence_lengths, index_to_endpoint_dict, known_adjacencies):
    """
    Split a scaffold path into two approximately equal halves.
    
    Parameters:
    -----------
    path : list
        Path indices representing a scaffold
    sequence_lengths : dict
        Sequence length information
    index_to_endpoint_dict : dict
        Mapping from indices to endpoint identifiers
    known_adjacencies : set
        Set of known adjacent endpoint pairs
        
    Returns:
    --------
    tuple
        (left_half, right_half)
    """
    path_length = len(path)
    assert path_length % 2 == 0

    ordered_components = list()

    for n in range(path_length // 2):
        endpoint_1 = index_to_endpoint_dict[path[2*n]]
        endpoint_2 = index_to_endpoint_dict[path[2*n+1]]

        comps_1, left_1, right_1 = extract_endpoint_components(endpoint_1)
        comps_2, left_2, right_2 = extract_endpoint_components(endpoint_2)

        # Determine orientation based on known adjacencies
        if tuple(sorted([left_1, left_2])) in known_adjacencies:
            ordered_components.extend(comps_1[::-1])
            ordered_components.extend(comps_2)
        elif tuple(sorted([right_1, right_2])) in known_adjacencies:
            ordered_components.extend(comps_1)
            ordered_components.extend(comps_2[::-1])
        elif tuple(sorted([left_1, right_2])) in known_adjacencies:
            ordered_components.extend(comps_1[::-1])
            ordered_components.extend(comps_2[::-1])
        else:
            assert tuple(sorted([right_1, left_2])) in known_adjacencies
            ordered_components.extend(comps_1)
            ordered_components.extend(comps_2)

    # Calculate midpoint for partitioning
    scaffold_tuple = tuple([index_to_endpoint_dict[i] for i in path])
    half_length = calculate_endpoint_length(scaffold_tuple, sequence_lengths) / 2

    accumulated = 0
    length_differences = list()

    for pos, component in enumerate(ordered_components):
        assert isinstance(component, str)
        accumulated += calculate_endpoint_length(component, sequence_lengths)
        length_differences.append((pos, abs(accumulated - half_length)))

    length_differences.sort(key=lambda x: x[1])
    split_position = length_differences[0][0] + 1

    left_half = tuple(ordered_components[:split_position])
    right_half = tuple(ordered_components[split_position:])

    # Record new adjacency
    adjacent_pair = tuple(sorted([left_half[-1], right_half[0]]))
    if adjacent_pair not in known_adjacencies:
        known_adjacencies.add(adjacent_pair)

    return left_half, right_half


def flatten_endpoint(endpoint):
    """
    Convert endpoint to list format regardless of input type.
    
    Parameters:
    -----------
    endpoint : str or tuple
        Endpoint identifier
        
    Returns:
    --------
    list
        List of component identifiers
    """
    if isinstance(endpoint, str):
        return [endpoint]
    else:
        assert isinstance(endpoint, tuple)
        return list(endpoint)


def update_graph_state(path_list, previous_contact_matrix, index_to_endpoint_dict,
                       endpoint_to_index_dict, flanking_dict, sequence_lengths,
                       known_adjacencies, flank_length):
    """
    Update graph state for next iteration.
    
    Parameters:
    -----------
    path_list : list
        List of paths from current iteration
    previous_contact_matrix : numpy.ndarray
        Contact matrix from previous iteration
    index_to_endpoint_dict : dict
        Current mapping from indices to endpoints
    endpoint_to_index_dict : dict
        Original mapping from endpoints to indices
    flanking_dict : dict
        Flanking region information
    sequence_lengths : dict
        Sequence length information
    known_adjacencies : set
        Set of known adjacent endpoint pairs
    flank_length : int
        Flanking region length
        
    Returns:
    --------
    tuple
        Updated (index_to_endpoint_dict, contact_dict, index_pairs, output_paths)
    """
    def add_to_output_path(endpoint):
        if isinstance(endpoint, tuple):
            output_paths[-1].extend(endpoint)
        else:
            assert isinstance(endpoint, str)
            output_paths[-1].append(endpoint)

    def compute_flanking_region(path_segment, direction):
        remaining_length = calculate_endpoint_length(path_segment, sequence_lengths)
        if remaining_length > flank_length:
            for idx, component in enumerate(path_segment[::direction]):
                comp_length = calculate_endpoint_length(component, sequence_lengths)
                if remaining_length - comp_length > flank_length:
                    remaining_length -= comp_length
                else:
                    break
            if idx == 0:
                flank_region = path_segment
            elif direction == 1:
                flank_region = path_segment[:-idx]
            else:
                flank_region = path_segment[idx:]
            flanking_dict[path_segment] = (flank_region, 
                                           calculate_endpoint_length(flank_region, sequence_lengths))

    index_pairs = list()
    new_index_to_endpoint = dict()
    output_paths = list()

    for path_idx, path in enumerate(path_list):
        idx_a, idx_b = 2 * path_idx, 2 * path_idx + 1
        index_pairs.append((idx_b, idx_a))
        output_paths.append([])

        if len(path) == 2:
            endpoint_left = index_to_endpoint_dict[path[0]]
            endpoint_right = index_to_endpoint_dict[path[1]]
            new_index_to_endpoint[idx_a] = endpoint_left
            new_index_to_endpoint[idx_b] = endpoint_right
            add_to_output_path(endpoint_left)
            add_to_output_path(endpoint_right)
        else:
            left_half, right_half = partition_scaffold_path(
                path, sequence_lengths, index_to_endpoint_dict, known_adjacencies)
            new_index_to_endpoint[idx_a] = left_half
            new_index_to_endpoint[idx_b] = right_half
            output_paths[-1].extend(left_half)
            output_paths[-1].extend(right_half)
            
            if flank_length:
                compute_flanking_region(left_half, -1)
                compute_flanking_region(right_half, 1)

    # Update contact dictionary
    contact_dict = defaultdict(int)

    for idx_i in new_index_to_endpoint:
        for idx_j in new_index_to_endpoint:
            if idx_i <= idx_j or (idx_i, idx_j) in index_pairs:
                continue

            endpoint_i = new_index_to_endpoint[idx_i]
            endpoint_j = new_index_to_endpoint[idx_j]

            if endpoint_i in flanking_dict:
                components_i = flatten_endpoint(flanking_dict[endpoint_i][0])
            else:
                components_i = flatten_endpoint(endpoint_i)
                
            if endpoint_j in flanking_dict:
                components_j = flatten_endpoint(flanking_dict[endpoint_j][0])
            else:
                components_j = flatten_endpoint(endpoint_j)

            for comp_i, comp_j in product(components_i, components_j):
                orig_i = endpoint_to_index_dict[comp_i]
                orig_j = endpoint_to_index_dict[comp_j]
                contact_count = previous_contact_matrix[orig_i, orig_j]
                if contact_count:
                    contact_dict[(idx_i, idx_j)] += contact_count

    return new_index_to_endpoint, contact_dict, index_pairs, output_paths


def write_tour_file(output_paths, prefix):
    """
    Write tour file for ALLHiC compatibility.
    
    Parameters:
    -----------
    output_paths : list
        List of output paths
    prefix : str
        Group prefix for output filename
    """
    tour_filename = '{}.tour'.format(prefix)
    with open(tour_filename, 'w') as tour_file:
        tour_file.write('>INIT\n')
        ordered_sequences = list()
        for path in output_paths:
            for endpoint in path[::2]:
                seq_name, orientation = endpoint.rsplit('_', 1)
                if orientation == 'H':
                    ordered_sequences.append(seq_name + '+')
                else:
                    ordered_sequences.append(seq_name + '-')
        tour_file.write('{}\n'.format(' '.join(ordered_sequences)))


def trim_shortest_path(index_pairs, contact_dict, density_matrix):
    """
    Remove the shortest path from current graph state.
    
    Parameters:
    -----------
    index_pairs : list
        List of index pairs
    contact_dict : dict
        Contact dictionary
    density_matrix : numpy.ndarray
        Density matrix
        
    Returns:
    --------
    numpy.ndarray
        Trimmed density matrix
    """
    idx_a, idx_b = index_pairs.pop(-1)

    for key in list(contact_dict.keys()):
        if idx_a in key or idx_b in key:
            contact_dict.pop(key)

    return density_matrix[:-2, :-2]


def execute_fast_ordering(args, sequence_lengths, group_data, prefix):
    """
    Perform fast ordering of sequences within a group.
    
    Parameters:
    -----------
    args : argparse.Namespace
        Command line arguments
    sequence_lengths : dict
        Sequence length information
    group_data : tuple
        Group-specific data (sequence_info, sequence_list, contact_dict, endpoint_index_map)
    prefix : str
        Group prefix
        
    Returns:
    --------
    tuple
        (output_paths, single_sequence_flag)
    """
    logger.info('[{}] Initiating fast ordering process...'.format(prefix))

    sequence_info, sequence_list, contact_dict, endpoint_index_map = group_data

    # Validate input data
    logger.info('[{}] Validating group file contents...'.format(prefix))
    for seq_name, seq_len in sequence_info:
        if seq_name not in sequence_lengths:
            logger.error('[{}] Sequence {} not found in FASTA file'.format(prefix, seq_name))
            raise RuntimeError('Group validation failed: [{}] Missing sequence {}'.format(prefix, seq_name))
        elif seq_len != sequence_lengths[seq_name]:
            logger.error('[{}] Length mismatch for sequence {}: group file reports {}, FASTA has {}'.format(
                prefix, seq_name, seq_len, sequence_lengths[seq_name]))
            raise RuntimeError('Group validation failed: [{}] Length mismatch for sequence {}'.format(prefix, seq_name))

    if len(sequence_info) == 1:
        logger.info('[{}] Single sequence group - ordering/orientation not required'.format(prefix))
        single_seq = sequence_info[0][0]
        return [[single_seq + '_H', single_seq + '_T']], True
    elif not sequence_info:
        logger.error('[{}] Empty group file detected'.format(prefix))
        raise RuntimeError('Group validation failed: [{}] No sequences found'.format(prefix))

    # Initialize data structures
    matrix_dimension = len(sequence_list) * 2
    index_to_endpoint = {idx: endpoint for endpoint, idx in endpoint_index_map.items()}

    index_pairs = list()
    output_paths = list()
    known_adjacencies = set()
    
    for seq in sequence_list:
        seq_h, seq_t = seq + '_H', seq + '_T'
        output_paths.append([seq_h, seq_t])
        index_pairs.append((endpoint_index_map[seq_h], endpoint_index_map[seq_t]))
        known_adjacencies.add(tuple(sorted([seq_h, seq_t])))

    current_paths = sequence_list.copy()
    flanking_regions = dict()
    removed_paths = list()
    skip_density = False
    iteration = 0

    logger.info('[{}] Beginning fast ordering iterations...'.format(prefix))

    while len(current_paths) != 1:
        iteration += 1

        if not skip_density:
            contact_matrix = sparse_matrix_converter(contact_dict, matrix_dimension)
            
            if iteration == 1:
                original_contact_matrix = contact_matrix

            density_matrix = compute_density_matrix(
                contact_matrix, matrix_dimension, index_to_endpoint, 
                sequence_lengths, flanking_regions, args.density_cal_method)

        confidence_matrix, max_density = calculate_edge_confidence(
            matrix_dimension, contact_dict, density_matrix)

        if max_density <= args.confidence_cutoff and len(current_paths) > 2:
            matrix_dimension -= 2
            density_matrix = trim_shortest_path(index_pairs, contact_dict, density_matrix)
            removed_paths.append(output_paths.pop(-1))
            current_paths.pop(-1)
            logger.debug('[{}] Iteration {}, MAXS {}'.format(prefix, iteration, max_density))
            logger.debug('[{}] Path {}, Length {}, removed'.format(
                prefix, len(output_paths) + 1, 
                calculate_endpoint_length(tuple(removed_paths[-1]), sequence_lengths)))
            skip_density = True
            continue
        elif max_density <= args.confidence_cutoff:
            assert len(current_paths) == 2
            break

        skip_density = False
        apply_confidence_filter(confidence_matrix, contact_dict, args.confidence_cutoff)

        # Build maximum spanning forest
        forest = nxtree.maximum_spanning_tree(Graph(confidence_matrix), algorithm='kruskal')

        current_paths.clear()

        for component_nodes in connected_components(forest):
            subtree = forest.subgraph(component_nodes)
            endpoints = []
            for node, degree in subtree.degree():
                if degree == 1:
                    endpoints.append(node)
            assert len(endpoints) == 2

            source, target = endpoints
            path = dict(shortest_path(subtree))[source][target]
            path_length = calculate_endpoint_length(
                tuple([index_to_endpoint[i] for i in path]), sequence_lengths)
            current_paths.append((path, path_length))

        # Sort by length descending
        current_paths = [path for path, length in 
                         sorted(current_paths, key=lambda x: x[1], reverse=True)]

        # Update for next iteration
        matrix_dimension = len(current_paths) * 2
        flank_region_size = args.flanking_region * 1000
        
        index_to_endpoint, contact_dict, index_pairs, output_paths = update_graph_state(
            current_paths, original_contact_matrix, index_to_endpoint,
            endpoint_index_map, flanking_regions, sequence_lengths,
            known_adjacencies, flank_region_size)

        logger.debug('[{}] Iteration {}, MAXS {}'.format(prefix, iteration, max_density))
        for idx, path in enumerate(output_paths, 1):
            logger.debug('[{}] Path {}, Length {}: {}'.format(
                prefix, idx, calculate_endpoint_length(tuple(path), sequence_lengths), 
                '->'.join(path)))

    output_paths.extend(removed_paths[::-1])
    return output_paths, False


def run_conhic_optimization(args, group_file, prefix, contact_map_path, allhic_path):
    """
    Execute ALLHiC optimization for a group.
    
    Parameters:
    -----------
    args : argparse.Namespace
        Command line arguments
    group_file : str
        Path to group file
    prefix : str
        Group prefix
    contact_map_path : str
        Path to contact map file
    allhic_path : str
        Path to ALLHiC executable
    """
    logger.info('[{}] Initiating ALLHiC optimization...'.format(prefix))

    cmd = [
        allhic_path, 'optimize', group_file, contact_map_path,
        '--mutapb', str(args.mutprob),
        '--ngen', str(args.ngen),
        '--npop', str(args.npop),
        '--seed', str(args.seed),
        '--resume'
    ]

    if args.skipGA:
        cmd.append('--skipGA')

    subprocess.run(
        cmd, check=True,
        stdout=open('allhic_logs/{}.allhic_stdout.log'.format(prefix), 'w'),
        stderr=open('allhic_logs/{}.allhic_stderr.log'.format(prefix), 'w')
    )


def compare_ordering_results(prefix, sequence_lengths):
    """
    Compare fast ordering results with ALLHiC optimization.
    
    Parameters:
    -----------
    prefix : str
        Group prefix
    sequence_lengths : dict
        Sequence length information
        
    Returns:
    --------
    bool
        True if fast ordering result is acceptable, False if ALLHiC result should be used
    """
    def longest_increasing_subsequence(values, length_dict, forward=True):
        if forward:
            filtered = [v for v in values if v > 0]
        else:
            filtered = [v for v in values if v < 0]

        if not filtered:
            return 0

        dp = [0] * len(filtered)
        predecessors = [None] * len(filtered)
        max_idx = 0

        for i in range(len(filtered)):
            dp[i] = length_dict[filtered[i]]
            for j in range(i):
                if filtered[i] > filtered[j] and dp[i] < dp[j] + length_dict[filtered[i]]:
                    dp[i] = dp[j] + length_dict[filtered[i]]
                    predecessors[i] = j
            if dp[i] >= dp[max_idx]:
                max_idx = i
        return dp[max_idx]

    fast_tour = '{}.tour.sav'.format(prefix)
    allhic_tour = '{}.tour'.format(prefix)

    fast_seqs, fast_oris = [], []
    for seq, ori in list(parse_tours([fast_tour], sequence_lengths)[0].values())[0]:
        fast_seqs.append(seq)
        fast_oris.append(ori)

    seq_lengths_list = [sequence_lengths[seq] for seq in fast_seqs]
    group_total = sum(seq_lengths_list)
    
    # Check if sequences are numerous enough to benefit from ALLHiC
    group_to_longest_ratio = group_total / max(seq_lengths_list)
    if group_to_longest_ratio > 50:
        logger.info('{}: Selecting ALLHiC optimization (group length / longest sequence = {})'.format(
            prefix, group_to_longest_ratio))
        return False

    allhic_seqs, allhic_oris = [], []
    for seq, ori in list(parse_tours([allhic_tour], sequence_lengths)[0].values())[0]:
        allhic_seqs.append(seq)
        allhic_oris.append(ori)

    best_lis_ratio = 0
    for rotation in range(len(fast_seqs) - 1):
        values, value_lengths = [], dict()
        for i, seq in enumerate(fast_seqs):
            j = allhic_seqs.index(seq)
            if fast_oris[i] == allhic_oris[j]:
                values.append(j + 1)
                value_lengths[j + 1] = sequence_lengths[seq]
            else:
                values.append(-j - 1)
                value_lengths[-j - 1] = sequence_lengths[seq]

        lis_forward = longest_increasing_subsequence(values, value_lengths, forward=True)
        lis_reverse = longest_increasing_subsequence(values, value_lengths, forward=False)
        max_lis = max(lis_forward, lis_reverse)
        lis_ratio = max_lis / group_total

        if lis_ratio >= 0.9:
            logger.info('{}: Selecting ALLHiC optimization (LIS length / group length = {})'.format(
                prefix, lis_ratio))
            return False
        else:
            if lis_ratio > best_lis_ratio:
                best_lis_ratio = lis_ratio
            
            # Rotate for next comparison
            fast_seqs = fast_seqs[1:] + [fast_seqs[0]]
            fast_oris = fast_oris[1:] + [fast_oris[0]]

    logger.info('{}: Selecting fast ordering (maximum LIS length / group length = {})'.format(
        prefix, best_lis_ratio))
    return True


def process_single_group(args, group_file, sequence_lengths, group_data, group_params, allhic_path):
    """
    Process a single group through the ordering pipeline.
    
    Parameters:
    -----------
    args : argparse.Namespace
        Command line arguments
    group_file : str
        Path to group file
    sequence_lengths : dict
        Sequence length information
    group_data : tuple
        Group-specific contact data
    group_params : tuple
        Group parameters (prefix, contact_map_path)
    allhic_path : str
        Path to ALLHiC executable
    """
    prefix, contact_map_path = group_params
    single_sequence = False
    use_fast_result = False

    output_paths, single_sequence = execute_fast_ordering(
        args, sequence_lengths, group_data, prefix)
    write_tour_file(output_paths, prefix)

    if not single_sequence:
        run_conhic_optimization(args, group_file, prefix, contact_map_path, allhic_path)
        use_fast_result = compare_ordering_results(prefix, sequence_lengths)

    final_tour_dir = 'final_tours'
    if use_fast_result:
        os.symlink('../{}.tour.sav'.format(prefix), 
                   os.path.join(final_tour_dir, '{}.tour'.format(prefix)))
    else:
        os.symlink('../{}.tour'.format(prefix), 
                   os.path.join(final_tour_dir, '{}.tour'.format(prefix)))


def check_multiprocessing_exceptions(result_list):
    """
    Check for exceptions in multiprocessing results.
    
    Parameters:
    -----------
    result_list : list
        List of AsyncResult objects
        
    Raises:
    -------
    RuntimeError
        If any exceptions were detected
    """
    error_count = 0
    for result in result_list:
        try:
            result.get()
        except Exception as e:
            error_count += 1
            logger.error(e)

    if error_count:
        raise RuntimeError('{} exception(s) detected - check logs above'.format(error_count))


def parse_command_line():
    """
    Parse command line arguments for ConHiC sort module.
    
    Returns:
    --------
    argparse.Namespace
        Parsed command line arguments
    """
    parser = argparse.ArgumentParser(prog='conhic sort')

    # Input file parameters
    input_section = parser.add_argument_group('Input File Parameters')
    input_section.add_argument(
        'fasta', help='Draft genome in FASTA format (use corrected assembly from clustering step)')
    input_section.add_argument(
        'contact_pickle', help='Contact dictionary pickle file (.pkl) from clustering step')
    input_section.add_argument(
        'contact_dir', help='Directory containing split contact map files (from reassignment step)')
    input_section.add_argument(
        'groups', nargs='+', help='Group files from reassignment step (group*.txt files)')

    # Fast ordering parameters
    fast_section = parser.add_argument_group('Fast Ordering Parameters')
    fast_section.add_argument(
        '--flanking_region', type=int, default=0,
        help='Use only flanking regions (ends) of sequences during fast ordering (kbp), default: 0 (whole sequences)')
    fast_section.add_argument(
        '--density_cal_method', choices={'multiplication', 'sum', 'geometric_mean'}, 
        default='multiplication',
        help='Method for contact density calculation, default: %(default)s')
    fast_section.add_argument(
        '--confidence_cutoff_range', type=str, default='1.01,1.03,0.01',
        help='Range of confidence cutoff values to test (start,end,step), default: %(default)s')

    # ALLHiC optimization parameters
    allhic_section = parser.add_argument_group('ALLHiC Optimization Parameters')
    allhic_section.add_argument(
        '--skipGA', default=False, action='store_true',
        help='Skip genetic algorithm optimization in ALLHiC, default: %(default)s')
    allhic_section.add_argument(
        '--mutprob', type=float, default=0.2,
        help='Mutation probability for genetic algorithm, default: %(default)s')
    allhic_section.add_argument(
        '--ngen', type=int, default=5000,
        help='Number of generations for convergence, default: %(default)s')
    allhic_section.add_argument(
        '--npop', type=int, default=100,
        help='Population size, default: %(default)s')
    allhic_section.add_argument(
        '--seed', type=int, default=42,
        help='Random seed, default: %(default)s')

    # Performance parameters
    perf_section = parser.add_argument_group('Performance Parameters')
    perf_section.add_argument(
        '--processes', type=int, default=24,
        help='Number of parallel processes (≤ number of group files), default: %(default)s')

    # Logging parameters
    log_section = parser.add_argument_group('Logging Parameters')
    log_section.add_argument(
        '--verbose', default=False, action='store_true',
        help='Enable verbose logging, default: %(default)s')

    return parser.parse_args()


def execute_single_cutoff(args, cutoff_value, log_filename=None):
    """
    Execute ConHiC sorting with a specific confidence cutoff.
    
    Parameters:
    -----------
    args : argparse.Namespace
        Command line arguments
    cutoff_value : float
        Confidence cutoff value for this run
    log_filename : str, optional
        Path to log file
    """
    if log_filename:
        file_handler = logging.FileHandler(log_filename, 'w')
        formatter = logging.Formatter(
            fmt='%(asctime)s <%(filename)s> [%(funcName)s] %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    start_time = time.time()
    logger.info('Running with confidence_cutoff = {}'.format(cutoff_value))
    logger.info('Program started, ConHiC version: {} (update: {})'.format(__version__, __update_time__))
    logger.info('Python version: {}'.format(sys.version.replace('\n', '')))
    logger.info('Command: {}'.format(' '.join(sys.argv)))

    process_count = min(args.processes, len(args.groups))

    if process_count > 1:
        pool = Pool(process_count)

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    # Locate ALLHiC executable
    logger.info('Checking ALLHiC executable location...')
    script_dir = os.path.dirname(os.path.realpath(__file__))
    allhic_exec = os.path.join(script_dir, 'allhic')

    if os.path.exists(allhic_exec):
        logger.info('ALLHiC found at {}'.format(script_dir))
    else:
        logger.error('ALLHiC not found at {}'.format(allhic_exec))
        raise RuntimeError('ALLHiC executable not found')

    # Load input data
    sequence_lengths = load_sequence_data(args.fasta)

    with open(args.contact_pickle, 'rb') as f:
        logger.info('Loading contact dictionary from pickle file...')
        full_contact_dict = pickle.load(f)

    # Extract group-specific data
    group_contact_data = dict()
    group_parameters = dict()

    logger.info('Processing group files and contact maps...')

    for group_file in args.groups:
        seq_info, contact_map, prefix = process_group_input(
            group_file, args.contact_dir)
        group_parameters[group_file] = (prefix, contact_map)
        seq_list = [seq for seq, length in seq_info]
        sub_contacts, endpoint_map = extract_subgraph_contacts(seq_list, full_contact_dict)
        group_contact_data[group_file] = (seq_info, seq_list, sub_contacts, endpoint_map)

    del full_contact_dict
    gc.collect()

    os.mkdir('allhic_logs')

    os.mkdir('final_tours')

    # Execute processing (parallel or sequential)
    if process_count > 1:
        logger.info('Using multiprocessing mode (processes={})'.format(process_count))
        results = []

        for group_file in args.groups:
            group_data = group_contact_data[group_file]
            group_contact_data[group_file] = None
            results.append(pool.apply_async(
                process_single_group, 
                args=(args, group_file, sequence_lengths, group_data, 
                      group_parameters[group_file], allhic_exec)))

        pool.close()
        pool.join()
        check_multiprocessing_exceptions(results)

    else:
        for group_file in args.groups:
            group_data = group_contact_data[group_file]
            group_contact_data[group_file] = None
            process_single_group(args, group_file, sequence_lengths, group_data,
                                 group_parameters[group_file], allhic_exec)

    elapsed = time.time() - start_time
    logger.info('Program completed in {:.2f}s'.format(elapsed))

    # Clean up file handlers
    for handler in logger.handlers[:]:
        if isinstance(handler, logging.FileHandler):
            logger.removeHandler(handler)
            handler.close()


def run(args, log_filename=None):
    """
    Main execution function for ConHiC sorting module.
    
    Parameters:
    -----------
    args : argparse.Namespace
        Command line arguments
    log_filename : str, optional
        Path to log file
    """
    # Parse confidence cutoff range
    cutoff_parts = args.confidence_cutoff_range.split(',')
    if len(cutoff_parts) != 3:
        raise ValueError('confidence_cutoff_range must be in format start,end,step')

    start_cutoff = float(cutoff_parts[0])
    end_cutoff = float(cutoff_parts[1])
    step_size = float(cutoff_parts[2])

    cutoff_values = []
    current = start_cutoff
    while current <= end_cutoff + 0.0001:
        cutoff_values.append(round(current, 2))
        current += step_size

    logger.info('Testing confidence cutoff values: {}'.format(cutoff_values))

    original_dir = os.getcwd()

    for cutoff in cutoff_values:
        logger.info('=' * 80)
        logger.info('Starting run with confidence_cutoff = {}'.format(cutoff))
        logger.info('=' * 80)

        cutoff_dir = os.path.join(original_dir, str(cutoff))
        os.makedirs(cutoff_dir, exist_ok=True)

        os.chdir(cutoff_dir)
        args.confidence_cutoff = cutoff

        execute_single_cutoff(args, cutoff, log_file='ConHiC_sort.log')

        os.chdir(original_dir)

        # Clean up for next iteration
        if os.path.exists('allhic_logs'):
            shutil.rmtree('allhic_logs')
        if os.path.exists('final_tours'):
            shutil.rmtree('final_tours')

        for filename in os.listdir('.'):
            if filename.endswith('.tour') or filename.endswith('.tour.sav'):
                os.remove(filename)

        logger.info('Completed run with confidence_cutoff = {}'.format(cutoff))
        logger.info('Results saved to: {}'.format(cutoff_dir))
        logger.info('-' * 80)


def main():
    """
    Entry point for ConHiC sorting module.
    """
    args = parse_command_line()
    run(args, log_file='ConHiC_sort.log')


if __name__ == '__main__':
    main()