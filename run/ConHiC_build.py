#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ConHiC: Chromosome-level scaffolding using Hi-C data
Build module - Constructs final scaffolds from ordering results
"""

import sys
import argparse
import os
import logging
import time
from collections import OrderedDict

# Local imports
from ConHiC_cluster import parse_fasta
from _version import __version__, __update_time__

# Configure logging
logging.basicConfig(
    format='%(asctime)s <%(filename)s> [%(funcName)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def parse_tours(tour_files, sequence_dict):
    """
    Parse ordering files containing scaffold component information
    
    Parameters:
    -----------
    tour_files : list
        List of .tour files from the ordering step
    sequence_dict : dict
        Dictionary containing sequence information
    
    Returns:
    --------
    tuple
        (tour_dict, assembled_sequences_set)
    """
    logger.info('Parsing tour files...')
    
    assembled_sequences = set()
    tour_dict = OrderedDict()
    
    for tour_file in tour_files:
        base_name = os.path.basename(tour_file)
        group_name = os.path.splitext(base_name)[0].rsplit('_', 1)[0]
        tour_dict[group_name] = list()
        
        last_line = ''
        with open(tour_file) as f:
            for line in f:
                if line.strip():
                    last_line = line.strip()
        
        for element in last_line.split():
            contig_id = element[:-1]
            orientation = element[-1]
            
            if contig_id not in sequence_dict:
                raise RuntimeError(f'Contig {contig_id} not found in FASTA file')
            elif contig_id in assembled_sequences:
                raise RuntimeError(f'Duplicate contig entry: {contig_id}')
            else:
                assembled_sequences.add(contig_id)
            
            tour_dict[group_name].append((contig_id, orientation))
    
    return tour_dict, assembled_sequences


def parse_corrected_ctgs(correction_file):
    """
    Load list of contigs that underwent sequence correction
    """
    corrected_set = set()
    
    if correction_file and os.path.exists(correction_file):
        with open(correction_file) as f:
            for line in f:
                if line.strip():
                    corrected_set.add(line.rstrip())
    
    return corrected_set


def build_final_scaffolds(tour_dict, sequence_dict, assembled_sequences, 
                          corrected_set, config):
    """
    Build final scaffold sequences and generate AGP files
    """
    gap_size = 100  # Default gap size
    
    def assemble_scaffold_sequence(group_name):
        """Assemble scaffold sequence by joining oriented contigs"""
        sequence_parts = []
        for contig_id, orientation in tour_dict[group_name]:
            if orientation == '+':
                sequence_parts.append(sequence_dict[contig_id][0])
            else:
                # Create complement translation table
                base_forward = 'ATCGNatcgn'
                base_reverse = 'TAGCNtagcn'
                com_tab = str.maketrans(base_forward, base_reverse)
                sequence_parts.append(
                    sequence_dict[contig_id][0].translate(com_tab)[::-1]
                )
        return ('N' * gap_size).join(sequence_parts)
    
    def write_agp_entries(group_name, agp_handle, raw_agp_handle):
        """Write AGP format entries for a scaffold group"""
        accumulated_length = 0
        
        for idx, (contig_id, orientation) in enumerate(tour_dict[group_name], 1):
            contig_length = sequence_dict[contig_id][1]
            start_pos = accumulated_length + 1
            end_pos = accumulated_length + contig_length
            accumulated_length += contig_length
            
            # Write contig component entry
            agp_handle.write(
                f'{group_name}\t{start_pos}\t{end_pos}\t{idx}\tW\t'
                f'{contig_id}\t1\t{contig_length}\t{orientation}\n'
            )
            
            # Write raw AGP entry
            if contig_id in corrected_set and ':' in contig_id:
                original_contig, coord_range = contig_id.rsplit(':', 1)
                start_coord, end_coord = coord_range.split('-')
                raw_agp_handle.write(
                    f'{group_name}\t{start_pos}\t{end_pos}\t{idx}\tW\t'
                    f'{original_contig}\t{start_coord}\t{end_coord}\t{orientation}\n'
                )
            else:
                raw_agp_handle.write(
                    f'{group_name}\t{start_pos}\t{end_pos}\t{idx}\tW\t'
                    f'{contig_id}\t1\t{contig_length}\t{orientation}\n'
                )
            
            # Add gap entry between contigs
            if idx < len(tour_dict[group_name]):
                gap_start = accumulated_length + 1
                gap_end = accumulated_length + gap_size
                accumulated_length += gap_size
                
                agp_handle.write(
                    f'{group_name}\t{gap_start}\t{gap_end}\t{idx+1}\tU\t'
                    f'{gap_size}\tscaffold\tyes\tproximity_ligation\n'
                )
                raw_agp_handle.write(
                    f'{group_name}\t{gap_start}\t{gap_end}\t{idx+1}\tU\t'
                    f'{gap_size}\tscaffold\tyes\tproximity_ligation\n'
                )
    
    logger.info('Building final scaffolds...')
    
    # Sort scaffolds by length (descending)
    scaffold_order = []
    for group, components in tour_dict.items():
        total_length = sum(sequence_dict[contig][1] for contig, _ in components) + \
                      (len(components) - 1) * gap_size
        scaffold_order.append((group, total_length))
    scaffold_order.sort(key=lambda x: x[1], reverse=True)
    scaffold_order = [group for group, _ in scaffold_order]
    
    # Process unassembled contigs
    unplaced_contigs = []
    for contig_id in sequence_dict:
        if contig_id not in assembled_sequences:
            unplaced_contigs.append((contig_id, sequence_dict[contig_id][1]))
    unplaced_contigs.sort(key=lambda x: x[1], reverse=True)
    
    # Write output files
    prefix = config.prefix if hasattr(config, 'prefix') else 'scaffolds'
    with open(f'{prefix}.fa', 'w') as fasta_out, \
         open(f'{prefix}.agp', 'w') as agp_out, \
         open(f'{prefix}.raw.agp', 'w') as raw_agp_out:
        
        # Process assembled scaffolds
        for group in scaffold_order:
            scaffold_seq = assemble_scaffold_sequence(group)
            fasta_out.write(f'>{group}\n{scaffold_seq}\n')
            write_agp_entries(group, agp_out, raw_agp_out)
        
        # Process unplaced contigs
        for contig_id, contig_length in unplaced_contigs:
            fasta_out.write(f'>{contig_id}\n{sequence_dict[contig_id][0]}\n')
            agp_out.write(f'{contig_id}\t1\t{contig_length}\t1\tW\t'
                         f'{contig_id}\t1\t{contig_length}\t+\n')
            
            if contig_id in corrected_set and ':' in contig_id:
                original_contig, coord_range = contig_id.rsplit(':', 1)
                start_coord, end_coord = coord_range.split('-')
                raw_agp_out.write(f'{contig_id}\t1\t{contig_length}\t1\tW\t'
                                 f'{original_contig}\t{start_coord}\t{end_coord}\t+\n')
            else:
                raw_agp_out.write(f'{contig_id}\t1\t{contig_length}\t1\tW\t'
                                 f'{contig_id}\t1\t{contig_length}\t+\n')


def generate_juicebox_script(config):
    """
    Generate shell script for Juicebox visualization
    """
    raw_fasta_name = os.path.basename(config.reference_assembly)
    script_directory = os.path.dirname(os.path.abspath(__file__))
    utils_directory = os.path.join(script_directory, '../utils')
    juicer_tool = os.path.join(utils_directory, 'juicer')
    juicer_jar = os.path.join(utils_directory, 'juicer_tools.1.9.9_jcuda.0.8.jar')
    prefix = config.prefix if hasattr(config, 'prefix') else 'scaffolds'
    
    with open('juicebox.sh', 'w') as script_file:
        script_file.write('#!/bin/bash\n\n')
        
        if not os.path.exists(raw_fasta_name):
            script_file.write(f'ln -s {config.reference_assembly} .\n')
        
        script_file.write(f'samtools faidx {raw_fasta_name}\n')
        script_file.write(
            f'{juicer_tool} pre -a -q 1 -o out_JBAT {config.mapping_data} '
            f'{prefix}.raw.agp {raw_fasta_name}.fai >out_JBAT.log 2>&1\n'
        )
        script_file.write(
            '(java -jar -Xmx32G {} pre out_JBAT.txt out_JBAT.hic.part '
            '<(cat out_JBAT.log | grep PRE_C_SIZE '.format(juicer_jar)
        )
        script_file.write(
            "| awk '{print $2\" \"$3}')) && (mv out_JBAT.hic.part out_JBAT.hic)\n"
        )


def parse_arguments():
    """
    Parse command-line arguments - simplified version with only essential parameters
    """
    parser = argparse.ArgumentParser(
        description='ConHiC build - Construct final scaffolds from ordering results'
    )
    
    parser.add_argument(
        'fasta',
        help='Draft genome in FASTA format'
    )
    
    parser.add_argument(
        'raw_fasta',
        help='Original draft genome for Juicebox visualization'
    )
    
    parser.add_argument(
        'alignments',
        help='Filtered Hi-C alignments in BAM format'
    )
    
    parser.add_argument(
        'tours', nargs='+',
        help='Tour files from the ordering step'
    )
    
    parser.add_argument(
        '--prefix', default='scaffolds',
        help='Output file prefix (default: scaffolds)'
    )
    
    # Hidden parameters with defaults (for internal use)
    parser.add_argument('--corrected_ctgs', default=None, help=argparse.SUPPRESS)
    parser.add_argument('--Ns', type=int, default=100, help=argparse.SUPPRESS)
    parser.add_argument('--sort_by_input', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--keep_letter_case', action='store_true', help=argparse.SUPPRESS)
    
    return parser.parse_args()


def run(args, log_file=None):
    """
    Main execution function (kept original name)
    """
    # Configure file logging if requested
    if log_file:
        file_handler = logging.FileHandler(log_file, 'w')
        formatter = logging.Formatter(
            fmt='%(asctime)s <%(filename)s> [%(funcName)s] %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    
    start_time = time.time()
    logger.info(f'ConHiC version: {__version__} (updated: {__update_time__})')
    logger.info(f'Command: {" ".join(sys.argv)}')
    
    # Simple parameter check
    if args.fasta != args.raw_fasta and not args.corrected_ctgs:
        logger.warning('Different assembly files provided but --corrected_ctgs is missing')
    
    # Load sequence data
    sequence_data = parse_fasta(
        args.fasta, 
        keep_letter_case=getattr(args, 'keep_letter_case', False),
        logger=logger
    )
    
    # Load corrected contigs if available
    corrected_file = getattr(args, 'corrected_ctgs', None)
    corrected_contigs = parse_corrected_ctgs(corrected_file)
    
    # Parse tour files
    tour_data, assembled_ids = parse_tours(args.tours, sequence_data)
    
    # Build scaffolds
    build_final_scaffolds(tour_data, sequence_data, assembled_ids, corrected_contigs, args)
    
    # Generate visualization script
    generate_juicebox_script(args)
    
    elapsed_time = time.time() - start_time
    logger.info(f'Processing completed in {elapsed_time:.2f} seconds')


def main():
    """
    Main entry point
    """
    args = parse_arguments()
    run(args, 'ConHiC_build.log')


if __name__ == '__main__':
    main()