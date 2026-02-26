import os
import sys
import subprocess
import argparse
import logging
import re
import glob
import time
import numpy as np
import ConHiC_cluster
import ConHiC_reassign as reassign_module
import ConHiC_build as build_module
import comprehensive_grouping as grouping_module
from _version import __version__, __update_time__

# Configure logging system
logging.basicConfig(
    format='%(asctime)s <%(filename)s> [%(funcName)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# Path utilities
abspath = os.path.abspath
joinpath = os.path.join
splitpath = os.path.split


def parse_arguments():
    """Parse command line arguments for the ConHiC pipeline"""
    
    parser = argparse.ArgumentParser(
        prog='conhic',
        description='ConHiC: Chromosome-scale scaffolding using Hi-C data'
    )

    # ========== Input file specifications ==========
    input_group = parser.add_argument_group('Input data specifications')
    input_group.add_argument(
        'assembly', 
        help='Input assembly in FASTA format'
    )
    input_group.add_argument(
        'hic_alignments', 
        help='Hi-C read alignments (BAM/pairs format, unsorted)'
    )
    input_group.add_argument(
        'chromosome_count', 
        type=int, 
        help='Expected number of chromosomes'
    )
    input_group.add_argument(
        '--format', 
        choices={'bam', 'pairs', 'bgzipped_pairs', 'auto'}, 
        default='auto',
        help='Alignment file format (default: %(default)s)'
    )
    input_group.add_argument(
        '--phasing_gfa', 
        default=None,
        help='Phased assembly GFA files (comma-separated) for haplotype separation'
    )
    input_group.add_argument(
        '--long_reads', 
        default=None,
        help='Ultra-long read alignments in BAM format'
    )

    # ========== Pipeline control ==========
    pipe_group = parser.add_argument_group('Pipeline control')
    pipe_group.add_argument(
        '--enzyme_sites', 
        default='GATC',
        help='Restriction enzyme recognition sites (comma-separated for multiple enzymes)'
    )
    pipe_group.add_argument(
        '--run_stages', 
        default='1,2,3,4',
        help='Pipeline stages to execute: 1=cluster, 2=reassign, 3=order, 4=build'
    )
    pipe_group.add_argument(
        '--quick_mode', 
        action='store_true',
        help='Skip clustering and reassignment, perform rapid ordering only'
    )
    pipe_group.add_argument(
        '--normalize_links', 
        action='store_true',
        help='Normalize Hi-C contacts by total links per contig/group'
    )
    pipe_group.add_argument(
        '--output_dir', 
        default=None,
        help='Output directory (default: current directory)'
    )
    pipe_group.add_argument(
        '--contig_lists', 
        default=[], 
        help='Pre-defined contig grouping lists'
    )

    # ========== Assembly refinement ==========
    refine_group = parser.add_argument_group('Assembly refinement')
    refine_group.add_argument(
        '--refine_rounds', 
        type=int, 
        default=0,
        help='Maximum refinement iterations (0 to disable)'
    )
    refine_group.add_argument(
        '--refine_resolution', 
        type=int, 
        default=500,
        help='Resolution for breakpoint detection (bp)'
    )
    refine_group.add_argument(
        '--coverage_ratio', 
        type=float, 
        default=0.2,
        help='Coverage threshold ratio for breakpoint identification'
    )
    refine_group.add_argument(
        '--region_length_ratio', 
        type=float, 
        default=0.1,
        help='Minimum region length ratio for high-coverage regions'
    )
    refine_group.add_argument(
        '--min_region_length', 
        type=int, 
        default=5000,
        help='Absolute minimum region length (bp)'
    )

    # ========== Preprocessing filters ==========
    filter_group = parser.add_argument_group('Preprocessing filters')
    filter_group.add_argument(
        '--length_percentile', 
        type=int, 
        default=80,
        help='Keep contigs longer than this percentile'
    )
    filter_group.add_argument(
        '--min_enzyme_sites', 
        type=int, 
        default=5,
        help='Minimum restriction enzyme sites per contig'
    )
    filter_group.add_argument(
        '--density_lower_limit', 
        default='0.2X',
        help='Lower limit for contact density filtering'
    )
    filter_group.add_argument(
        '--density_upper_limit', 
        default='1.9X',
        help='Upper limit for contact density filtering'
    )
    filter_group.add_argument(
        '--depth_upper_limit', 
        default='1.5X',
        help='Upper limit for read depth filtering'
    )
    filter_group.add_argument(
        '--neighbor_count', 
        type=int, 
        default=10,
        help='Number of nearest neighbors for rank-sum calculation'
    )
    filter_group.add_argument(
        '--rank_sum_cutoff', 
        type=int, 
        default=0,
        help='Hard cutoff for rank-sum values'
    )
    filter_group.add_argument(
        '--rank_sum_upper_limit', 
        default='1.5X',
        help='Upper limit for rank-sum filtering'
    )
    filter_group.add_argument(
        '--ploidy', 
        type=int, 
        default=0,
        help='Remove inter-allelic Hi-C contacts (ploidy >= 2)'
    )
    filter_group.add_argument(
        '--concordance_ratio', 
        type=float, 
        default=0.2,
        help='Concordance threshold for allelic pair identification'
    )
    filter_group.add_argument(
        '--window_count', 
        type=int, 
        default=50,
        help='Number of windows for concordance calculation'
    )
    filter_group.add_argument(
        '--filter_concentrated', 
        action='store_true',
        help='Remove locally concentrated Hi-C contacts'
    )
    filter_group.add_argument(
        '--max_read_pairs', 
        type=int, 
        default=200,
        help='Maximum read pairs for allelic analysis'
    )
    filter_group.add_argument(
        '--min_read_pairs', 
        type=int, 
        default=20,
        help='Minimum read pairs for allelic analysis'
    )
    filter_group.add_argument(
        '--phasing_weight', 
        type=float, 
        default=1.0,
        help='Weight for phasing information (0-1)'
    )

    # ========== Long read processing ==========
    longread_group = parser.add_argument_group('Long read processing')
    longread_group.add_argument(
        '--min_longread_mapq', 
        type=int, 
        default=30,
        help='Minimum MAPQ for long read alignments'
    )
    longread_group.add_argument(
        '--min_longread_length', 
        type=int, 
        default=10000,
        help='Minimum alignment length for long reads'
    )
    longread_group.add_argument(
        '--max_end_distance', 
        type=int, 
        default=100,
        help='Maximum distance to contig ends'
    )
    longread_group.add_argument(
        '--max_overlap_ratio', 
        type=float, 
        default=0.5,
        help='Maximum overlap ratio between alignments'
    )
    longread_group.add_argument(
        '--max_gap', 
        type=int, 
        default=10000,
        help='Maximum gap between alignments'
    )
    longread_group.add_argument(
        '--min_longread_support', 
        type=int, 
        default=2,
        help='Minimum long reads supporting a contig pair'
    )

    # ========== Clustering parameters ==========
    cluster_group = parser.add_argument_group('Clustering parameters')
    cluster_group.add_argument(
        '--window_size', 
        type=int, 
        default=-1,
        help='Window size for matrix construction (Kbp)'
    )
    cluster_group.add_argument(
        '--end_region', 
        type=int, 
        default=500,
        help='Use only contig ends for matrix construction (Kbp)'
    )
    cluster_group.add_argument(
        '--expansion', 
        type=int, 
        default=2,
        help='Expansion parameter for MCL'
    )
    cluster_group.add_argument(
        '--min_inflation', 
        type=float, 
        default=1.1,
        help='Minimum inflation for MCL'
    )
    cluster_group.add_argument(
        '--max_inflation', 
        type=float, 
        default=1.2,
        help='Maximum inflation for MCL'
    )
    cluster_group.add_argument(
        '--inflation_step', 
        type=float, 
        default=0.1,
        help='Inflation step size'
    )
    cluster_group.add_argument(
        '--max_iterations', 
        type=int, 
        default=200,
        help='Maximum MCL iterations'
    )
    cluster_group.add_argument(
        '--pruning_threshold', 
        type=float, 
        default=0.0001,
        help='Pruning threshold for MCL'
    )
    cluster_group.add_argument(
        '--skip_clustering', 
        action='store_true',
        help='Skip MCL clustering'
    )

    # ========== Reassignment parameters ==========
    reassign_group = parser.add_argument_group('Reassignment parameters')
    reassign_group.add_argument(
        '--min_group_size', 
        type=float, 
        default=5,
        help='Minimum group size (Mbp)'
    )
    reassign_group.add_argument(
        '--max_reassign_length', 
        type=float, 
        default=10000,
        help='Maximum contig length for reassignment (Kbp)'
    )
    reassign_group.add_argument(
        '--min_reassign_enzyme', 
        type=int, 
        default=25,
        help='Minimum enzyme sites for reassignment'
    )
    reassign_group.add_argument(
        '--min_reassign_links', 
        type=int, 
        default=25,
        help='Minimum Hi-C links for reassignment'
    )
    reassign_group.add_argument(
        '--min_reassign_density', 
        type=float, 
        default=0.0001,
        help='Minimum contact density for reassignment'
    )
    reassign_group.add_argument(
        '--density_ratio', 
        type=float, 
        default=4,
        help='Minimum density ratio for reassignment'
    )
    reassign_group.add_argument(
        '--ambiguity_threshold', 
        type=float, 
        default=0.6,
        help='Threshold for ambiguous contigs'
    )
    reassign_group.add_argument(
        '--reassign_rounds', 
        type=int, 
        default=5,
        help='Maximum reassignment rounds'
    )
    reassign_group.add_argument(
        '--skip_rescue', 
        action='store_true',
        help='Skip additional rescue round'
    )

    # ========== Ordering parameters ==========
    order_group = parser.add_argument_group('Ordering parameters')
    order_group.add_argument(
        '--skip_fast_order', 
        action='store_true',
        help='Skip fast ordering step'
    )
    order_group.add_argument(
        '--order_end_region', 
        type=int, 
        default=0,
        help='Use contig ends for ordering (Kbp)'
    )
    order_group.add_argument(
        '--density_method', 
        choices={'multiplication', 'sum', 'geometric_mean'}, 
        default='multiplication',
        help='Method for contact density calculation'
    )
    order_group.add_argument(
        '--confidence_threshold', 
        type=float, 
        default=1,
        help='Confidence threshold for ordering'
    )

    # ========== ALLHiC optimization ==========
    allhic_group = parser.add_argument_group('ALLHiC optimization')
    allhic_group.add_argument(
        '--skip_allhic', 
        action='store_true',
        help='Skip ALLHiC optimization'
    )
    allhic_group.add_argument(
        '--skip_ga', 
        action='store_true',
        help='Skip genetic algorithm'
    )
    allhic_group.add_argument(
        '--mutation_prob', 
        type=float, 
        default=0.2,
        help='Mutation probability in GA'
    )
    allhic_group.add_argument(
        '--generations', 
        type=int, 
        default=5000,
        help='Number of GA generations'
    )
    allhic_group.add_argument(
        '--population_size', 
        type=int, 
        default=100,
        help='Population size in GA'
    )
    allhic_group.add_argument(
        '--random_seed', 
        type=int, 
        default=42,
        help='Random seed for reproducibility'
    )

    # ========== Output parameters ==========
    output_group = parser.add_argument_group('Output parameters')
    output_group.add_argument(
        '--gap_fill', 
        type=int, 
        default=100,
        help='Number of Ns for gap filling'
    )
    output_group.add_argument(
        '--sort_by_input', 
        action='store_true',
        help='Sort scaffolds by input order'
    )
    output_group.add_argument(
        '--preserve_case', 
        action='store_true',
        help='Preserve original base case'
    )
    output_group.add_argument(
        '--output_prefix', 
        default='scaffolds',
        help='Prefix for output files'
    )

    # ========== Performance parameters ==========
    perf_group = parser.add_argument_group('Performance parameters')
    perf_group.add_argument(
        '--threads', 
        type=int, 
        default=8,
        help='Number of threads for I/O operations'
    )
    perf_group.add_argument(
        '--use_dense', 
        action='store_true',
        help='Use dense matrices instead of sparse'
    )
    perf_group.add_argument(
        '--processes', 
        type=int, 
        default=8,
        help='Number of parallel processes'
    )

    # ========== Logging ==========
    log_group = parser.add_argument_group('Logging')
    log_group.add_argument(
        '--verbose', 
        action='store_true',
        help='Enable verbose logging'
    )

    return parser.parse_args()


def expand_enzyme_sites(sites):
    """Expand degenerate restriction enzyme sites"""
    
    expanded_sites = []
    
    for site in sites:
        if 'N' in site:
            expanded_sites.append(site.replace('N', 'A', 1))
            expanded_sites.append(site.replace('N', 'T', 1))
            expanded_sites.append(site.replace('N', 'C', 1))
            expanded_sites.append(site.replace('N', 'G', 1))
        else:
            expanded_sites.append(site)
    
    if 'N' not in ''.join(expanded_sites):
        return expanded_sites
    else:
        return expand_enzyme_sites(expanded_sites)


def count_enzyme_sites(sequence, enzymes):
    """Count restriction enzyme sites in a sequence"""
    
    sites = [site.strip().upper() for site in enzymes.split(',') if site.strip()]
    expanded_sites = expand_enzyme_sites(sites)
    
    site_count = 0
    for site in expanded_sites:
        site_count += sequence.count(site)
    
    return site_count


def parse_assembly(assembly_file, enzymes='GATC', preserve_case=False):
    """Parse assembly FASTA and extract sequence information"""
    
    logger.info('Reading assembly file...')
    
    assembly_dict = {}
    with open(assembly_file) as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith('>'):
                contig = line.split()[0][1:]
                assembly_dict[contig] = []
            else:
                if preserve_case:
                    assembly_dict[contig].append(line.strip())
                else:
                    assembly_dict[contig].append(line.strip().upper())
    
    for contig, seq_parts in assembly_dict.items():
        sequence = ''.join(seq_parts)
        # Add pseudo-count of 1 to avoid division by zero
        enzyme_count = count_enzyme_sites(sequence, enzymes) + 1
        assembly_dict[contig] = [sequence, len(sequence), enzyme_count]
    
    return assembly_dict


def stage1_clustering(config):
    """Execute first pipeline stage: contig clustering"""
    
    logger.info('Stage 1: Clustering contigs using Hi-C contacts...')
    
    work_dir = '01.clustering'
    log_file = 'clustering.log'
    random_rounds = 10
    
    if not os.path.exists(work_dir):
        os.mkdir(work_dir)
    
    os.chdir(work_dir)
    
    # Run clustering algorithm
    cluster_results, window_sizes = ConHiC_cluster.run(
        config, 
        log_file=log_file,
        random_rounds=random_rounds
    )
    
    logger.info(f'Clustering results: {cluster_results}')
    
    # Update config for subsequent stages
    if config.refine_rounds:
        config.assembly = abspath('corrected_assembly.fa')
        config.corrected_contigs = abspath('corrected_contigs.txt')
        if config.quick_mode and config.phasing_gfa and len(config.phasing_gfa.split(',')) >= 2:
            config.phasing_gfa = ','.join([
                joinpath(splitpath(gfa)[0], work_dir, 'corrected_' + splitpath(gfa)[1]) 
                for gfa in config.phasing_gfa.split(',')
            ])
    else:
        config.corrected_contigs = None
    
    config.corrected_contigs = None
    os.chdir('..')
    
    # Prepare clustering outputs for multiple window sizes
    config.links = config.hic_links = config.clm_files = []
    config.cluster_sets = {}
    config.cluster_count = config.chromosome_count
    
    # Analyze contig length distribution
    assembly_info = parse_assembly(config.assembly, enzymes=config.enzyme_sites)
    logger.info(f'{len(assembly_info)} contigs in assembly')
    
    # Calculate window sizes based on contig lengths
    lengths = [item[1] for item in assembly_info.values()]
    mean_length = np.mean(lengths)
    
    window_sizes = [
        mean_length, mean_length*2, mean_length*1.5, mean_length/2,
        mean_length/4, mean_length/6, mean_length/8, mean_length/10
    ]
    window_sizes = [int(val) for val in window_sizes]
    config.window_size_list = window_sizes
    
    # Store paths for each window size
    config.links = joinpath(work_dir, str(int(mean_length)), 'full_contacts.pkl')
    config.hic_links = joinpath(work_dir, str(int(mean_length)), 'filtered_contacts.pkl')
    config.clm_files = joinpath(work_dir, str(int(mean_length)), 'paired_contacts.clm')
    
    for win_size in window_sizes:
        if win_size > int(mean_length/6):
            win_dir = joinpath(work_dir, str(win_size))
            inflation = cluster_results[str(win_size)]
            config.cluster_sets[str(win_size)] = joinpath(
                win_dir, f'inflation_{inflation}',
                f'mcl_inflation_{inflation}.clusters.txt'
            )
        else:
            config.cluster_sets[str(win_size)] = []
            for i in range(random_rounds):
                win_dir = joinpath(work_dir, str(win_size), str(i))
                inflation = cluster_results[str(win_size)]
                config.cluster_sets[str(win_size)].append(joinpath(
                    win_dir, f'inflation_{inflation}',
                    f'mcl_inflation_{inflation}.clusters.txt'
                ))


def stage2_reassignment(config):
    """Execute second pipeline stage: contig reassignment and rescue"""
    
    logger.info('Stage 2: Reassigning and rescuing contigs...')
    
    work_dir = '02.reassignment'
    log_file = 'reassignment.log'
    random_rounds = 10
    multi_scale_results = []
    
    if not os.path.exists(work_dir):
        os.mkdir(work_dir)
    
    os.chdir(work_dir)
    
    cluster_sets = config.cluster_sets
    
    for win_size in config.window_size_list:
        win_subdir = str(win_size)
        if not os.path.exists(win_subdir):
            os.mkdir(win_subdir)
        else:
            result_file = joinpath(work_dir, win_subdir, 'final_groups/final_clusters.txt')
            multi_scale_results.append(result_file)
            continue
        
        os.chdir(win_subdir)
        
        if win_size <= int(config.window_size_list[0]/6):
            cluster_sets[str(win_size)] = grouping_module.run(
                cluster_sets[str(win_size)]
            )
        
        reassign_module.run(
            config, 
            log_file=log_file,
            random_rounds=random_rounds,
            cluster_file=cluster_sets[str(win_size)]
        )
        
        os.chdir('..')
        result_file = joinpath(work_dir, win_subdir, 'final_groups/final_clusters.txt')
        multi_scale_results.append(result_file)
    
    # Combine results from all scales
    final_clusters = grouping_module.run(multi_scale_results)
    reassign_module.run(
        config, 
        log_file=log_file,
        random_rounds=random_rounds,
        cluster_file=final_clusters
    )
    
    os.chdir('..')
    
    # Update config for next stage
    config.clm_dir = abspath(joinpath(work_dir, 'split_clms'))
    config.group_files = [
        f for f in glob.glob(abspath(joinpath(work_dir, 'final_groups/group*.txt')))
    ]


def stage3_ordering(config):
    """Execute third pipeline stage: contig ordering and orientation"""
    
    logger.info('Stage 3: Ordering and orienting contigs within groups...')
    
    work_dir = '03.ordering'
    
    os.mkdir(work_dir)
    os.chdir(work_dir)
    
    # Call ordering script via subprocess for better memory management
    script_path = joinpath(
        os.path.dirname(os.path.realpath(__file__)), 
        'ConHiC_order.py'
    )
    
    command = [script_path, config.assembly, config.hic_links, config.clm_dir]
    command.extend(config.group_files)
    
    if config.quick_mode:
        command.append('--quick_mode')
    else:
        if config.skip_fast_order:
            command.append('--skip_fast_order')
        if config.skip_allhic:
            command.append('--skip_allhic')
        if config.skip_ga:
            command.append('--skip_ga')
        
        command.extend([
            '--mutprob', str(config.mutation_prob),
            '--ngen', str(config.generations),
            '--npop', str(config.population_size),
            '--seed', str(config.random_seed)
        ])
        command.extend(['--processes', str(config.processes)])
    
    command.extend(['--end_region', str(config.order_end_region)])
    command.extend(['--density_method', config.density_method])
    command.extend(['--confidence', str(config.confidence_threshold)])
    
    if config.verbose:
        command.append('--verbose')
    
    subprocess.run(command, check=True)
    
    os.chdir('..')
    
    # Update config for final stage
    config.tour_files = [
        f for f in glob.glob(abspath(joinpath(work_dir, 'final_tours', '*.tour')))
    ]


def stage4_build(config):
    """Execute final pipeline stage: scaffold construction"""
    
    logger.info('Stage 4: Building final chromosome-scale scaffolds...')
    
    work_dir = '04.build'
    log_file = 'build.log'
    
    os.mkdir(work_dir)
    os.chdir(work_dir)
    
    # Convert pairs format to BED if needed
    if config.hic_alignments.endswith('.pairs') or config.hic_alignments.endswith('.pairs.gz'):
        config.hic_alignments = '../01.clustering/alignments.bed'
    
    build_module.run(config, log_file=log_file)
    
    os.chdir('..')


def main():
    """Main execution function"""
    
    # Parse command line arguments
    config = parse_arguments()
    
    start_time = time.time()
    logger.info(f'ConHiC pipeline started (version: {__version__}, updated: {__update_time__})')
    logger.info(f'Python version: {sys.version.replace(chr(10), "")}')
    logger.info(f'Command: {" ".join(sys.argv)}')
    
    # Validate pipeline stages
    stages = {int(step) for step in config.run_stages.split(',')}
    valid_stages = [{1}, {1, 2}, {1, 2, 3}, {1, 2, 3, 4}]
    
    if stages not in valid_stages:
        raise ValueError(f'Invalid stage combination: {config.run_stages}')
    
    # Normalize paths
    config.assembly = abspath(config.assembly)
    config.original_assembly = config.assembly
    config.hic_alignments = abspath(config.hic_alignments)
    
    if config.long_reads:
        config.long_reads = abspath(config.long_reads)
    
    if config.phasing_gfa:
        gfa_files = [abspath(gfa) for gfa in config.phasing_gfa.split(',')]
        config.phasing_gfa = ','.join(gfa_files)
    
    # Setup output directory
    if config.output_dir:
        try:
            os.mkdir(config.output_dir)
        except FileExistsError:
            logger.warning(f'Directory already exists: {config.output_dir}')
        os.chdir(config.output_dir)
    
    # Execute pipeline stages
    stage1_clustering(config)
    
    if 2 in stages:
        stage2_reassignment(config)
    
    if 3 in stages:
        stage3_ordering(config)
    
    if 4 in stages:
        stage4_build(config)
    
    elapsed_time = time.time() - start_time
    logger.info(f'ConHiC pipeline completed in {elapsed_time:.2f}s')


if __name__ == '__main__':
    main()