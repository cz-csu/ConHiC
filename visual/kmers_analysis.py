#!/usr/bin/env python3
"""
Kmer-based genome assembly annotation tool
Annotates assembly contigs based on reference genome kmer composition
"""

import sys
import argparse
import logging
import collections
from typing import Dict, List, Tuple, Set, DefaultDict, Optional
from pathlib import Path

# Constants
COMPLEMENT_TABLE = str.maketrans('ATCG', 'TAGC')
CHROMOSOME_PREFIX = 'chr'
SHARED_SUFFIX = '_shared'
UNRELIABLE_CATEGORIES = {'nonspecific', 'unknown'}
VALID_HAP_TYPES = {f'hap{i}' for i in range(1, 5)}

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


class KmerAnnotationError(Exception):
    """Custom exception for kmer annotation errors"""
    pass


def read_fasta(file_path: str) -> Dict[str, str]:
    """
    Read FASTA file and return dictionary of sequences
    
    Args:
        file_path: Path to FASTA file
        
    Returns:
        Dictionary with sequence IDs as keys and sequences as values
        
    Raises:
        FileNotFoundError: If file doesn't exist
        ValueError: If FASTA format is invalid
    """
    if not Path(file_path).exists():
        raise FileNotFoundError(f"FASTA file not found: {file_path}")
    
    sequences: Dict[str, str] = {}
    current_id: Optional[str] = None
    current_seq: List[str] = []
    
    try:
        with open(file_path, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line:
                    continue
                    
                if line.startswith('>'):
                    # Save previous sequence
                    if current_id and current_seq:
                        sequences[current_id] = ''.join(current_seq).upper()
                        current_seq = []
                    
                    # Extract sequence ID
                    current_id = line.split()[0][1:]
                    if not current_id:
                        raise ValueError(f"Empty sequence ID at line {line_num}")
                else:
                    if current_id is None:
                        raise ValueError(f"Sequence without header at line {line_num}")
                    current_seq.append(line)
        
        # Save last sequence
        if current_id and current_seq:
            sequences[current_id] = ''.join(current_seq).upper()
            
    except IOError as e:
        raise IOError(f"Error reading file {file_path}: {e}")
    
    if not sequences:
        raise ValueError(f"No valid sequences found in {file_path}")
    
    logger.info(f"Loaded {len(sequences)} sequences from {file_path}")
    return sequences


def get_canonical_kmer(kmer: str) -> str:
    """
    Convert kmer to canonical form (lexicographically smaller of kmer and its reverse complement)
    
    Args:
        kmer: DNA kmer string
        
    Returns:
        Canonical kmer
    """
    if 'N' in kmer:
        return kmer
    reverse_comp = kmer.translate(COMPLEMENT_TABLE)[::-1]
    return kmer if kmer < reverse_comp else reverse_comp


def extract_kmers_from_sequence(seq: str, k: int) -> List[str]:
    """
    Extract all canonical kmers from a sequence
    
    Args:
        seq: DNA sequence
        k: Kmer size
        
    Returns:
        List of canonical kmers
    """
    kmers = []
    seq_len = len(seq)
    
    for i in range(seq_len - k + 1):
        kmer_forward = seq[i:i + k]
        if 'N' in kmer_forward:
            continue
        kmers.append(get_canonical_kmer(kmer_forward))
    
    return kmers


def build_kmer_index(sequences: Dict[str, str], k: int) -> Tuple[List[List[str]], List[str]]:
    """
    Build kmer index for all sequences
    
    Args:
        sequences: Dictionary of sequence IDs and sequences
        k: Kmer size
        
    Returns:
        Tuple of (list of kmer lists per sequence, list of sequence IDs)
    """
    kmer_lists = []
    seq_ids = []
    
    for seq_id, seq in sequences.items():
        logger.info(f"Extracting kmers from {seq_id}")
        kmers = extract_kmers_from_sequence(seq, k)
        kmer_lists.append(kmers)
        seq_ids.append(seq_id)
        logger.debug(f"Extracted {len(kmers)} kmers from {seq_id}")
    
    return kmer_lists, seq_ids


def annotate_kmers(
    ref_sequences: Dict[str, str],
    k: int,
    query_kmer_lists: List[List[str]]
) -> Dict[str, DefaultDict[str, int]]:
    """
    Annotate query kmers based on reference sequences
    
    Args:
        ref_sequences: Reference genome sequences
        k: Kmer size
        query_kmer_lists: List of kmer lists from query sequences
        
    Returns:
        Dictionary mapping kmers to source counts
    """
    # Build set of unique query kmers
    query_kmers = set()
    for kmers in query_kmer_lists:
        query_kmers.update(kmers)
    
    logger.info(f"Total unique query kmers: {len(query_kmers)}")
    
    # Annotate kmers from reference
    kmer_annotation: Dict[str, DefaultDict[str, int]] = {}
    
    for ref_id, ref_seq in ref_sequences.items():
        logger.info(f"Annotating kmers from reference {ref_id}")
        ref_kmers = extract_kmers_from_sequence(ref_seq, k)
        
        for kmer in ref_kmers:
            if kmer in query_kmers:
                if kmer not in kmer_annotation:
                    kmer_annotation[kmer] = collections.defaultdict(int)
                kmer_annotation[kmer][ref_id] += 1
    
    logger.info(f"Annotated {len(kmer_annotation)} kmers")
    return kmer_annotation


def classify_bin_kmers(
    kmers: List[str],
    kmer_annotation: Dict[str, DefaultDict[str, int]],
    estimated_chromosome: str
) -> Tuple[DefaultDict[str, int], float]:
    """
    Classify kmers in a genomic bin
    
    Args:
        kmers: List of kmers in the bin
        kmer_annotation: Kmer annotation dictionary
        estimated_chromosome: Estimated chromosome for this contig
        
    Returns:
        Tuple of (classification statistics, alpha value)
    """
    stats: DefaultDict[str, int] = collections.defaultdict(int)
    chr_specific = 0
    
    for kmer in kmers:
        if kmer not in kmer_annotation:
            stats['unreliable'] += 1
            continue
        
        sources = kmer_annotation[kmer]
        
        # Single source
        if len(sources) == 1:
            source = list(sources.keys())[0]
            if source.startswith(CHROMOSOME_PREFIX):
                chr_name, hap = source.split('_')
                hap = f'hap{hap}'
                if chr_name == estimated_chromosome:
                    stats[hap] += 1
                    chr_specific += 1
                else:
                    stats['other_chrom'] += 1
            else:
                stats['unreliable'] += 1
        
        # Multiple sources
        else:
            source_list = list(sources.keys())
            chr_names = {s.split('_')[0] for s in source_list if s.startswith(CHROMOSOME_PREFIX)}
            
            if len(chr_names) == 1:
                # Shared between haplotypes of same chromosome
                chr_name = chr_names.pop()
                if chr_name == estimated_chromosome:
                    stats['shared'] += 1
                    chr_specific += 1
                else:
                    stats['other_chrom'] += 1
            else:
                stats['unreliable'] += 1
    
    # Calculate alpha (confidence score)
    alpha = 1.0
    if chr_specific > 0:
        # Find maximum haplotype count
        max_hap = 0
        for cat in VALID_HAP_TYPES:
            if cat in stats and stats[cat] > max_hap:
                max_hap = stats[cat]
        if max_hap > 0:
            alpha = max_hap / chr_specific
    
    return stats, alpha


def determine_primary_source(stats: DefaultDict[str, int]) -> str:
    """
    Determine the primary source category for a bin
    
    Args:
        stats: Classification statistics
        
    Returns:
        Primary source category
    """
    if not stats:
        return 'unknown'
    
    # Sort categories by count
    sorted_cats = sorted(stats.items(), key=lambda x: x[1])
    primary = sorted_cats[-1][0]
    
    # For shared regions, check if a specific haplotype dominates
    if primary == 'shared':
        hap_counts = {cat: count for cat, count in stats.items() 
                     if cat in VALID_HAP_TYPES}
        if hap_counts:
            primary = max(hap_counts.items(), key=lambda x: x[1])[0]
    
    return primary


def process_assembly(
    assembly_seq: str,
    assembly_id: str,
    kmer_annotation: Dict[str, DefaultDict[str, int]],
    k: int,
    bin_size: int
) -> List[Dict]:
    """
    Process a single assembly sequence
    
    Args:
        assembly_seq: Assembly sequence
        assembly_id: Assembly sequence ID
        kmer_annotation: Kmer annotation dictionary
        k: Kmer size
        bin_size: Bin size for segmentation
        
    Returns:
        List of bin results
    """
    # Extract kmers from assembly
    assembly_kmers = extract_kmers_from_sequence(assembly_seq, k)
    
    # Estimate chromosome for this contig
    chr_counts: DefaultDict[str, int] = collections.defaultdict(int)
    for kmer in assembly_kmers:
        if kmer in kmer_annotation:
            for source in kmer_annotation[kmer]:
                if source.startswith(CHROMOSOME_PREFIX):
                    chr_name = source.split('_')[0]
                    chr_counts[chr_name] += 1
    
    estimated_chromosome = max(chr_counts.items(), key=lambda x: x[1])[0] if chr_counts else 'unknown'
    logger.info(f"Estimated chromosome for {assembly_id}: {estimated_chromosome}")
    
    # Process bins
    bin_results = []
    num_bins = (len(assembly_kmers) + bin_size - 1) // bin_size
    
    for bin_idx in range(num_bins):
        start = bin_idx * bin_size
        end = min(start + bin_size, len(assembly_kmers))
        bin_kmers = assembly_kmers[start:end]
        
        if not bin_kmers:
            continue
        
        stats, alpha = classify_bin_kmers(bin_kmers, kmer_annotation, estimated_chromosome)
        primary_source = determine_primary_source(stats)
        
        bin_results.append({
            'contig': assembly_id,
            'start': start * k + 1,  # Convert to genomic coordinates
            'end': end * k + k - 1,
            'primary_source': primary_source,
            'alpha': alpha,
            'stats': dict(stats)
        })
    
    return bin_results


def write_output(results: List[Dict], output_prefix: str):
    """
    Write results to output file
    
    Args:
        results: List of bin results
        output_prefix: Prefix for output filename
    """
    output_file = f"{output_prefix}_annotation.txt"
    
    with open(output_file, 'w') as f:
        # Write header
        f.write("#contig\tstart\tend\tprimary_source\talpha\tstatistics\n")
        
        for result in results:
            f.write(f"{result['contig']}\t{result['start']}\t{result['end']}\t"
                   f"{result['primary_source']}\t{result['alpha']:.3f}\t"
                   f"{result['stats']}\n")
    
    logger.info(f"Results written to {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Annotate genome assembly using kmer-based reference comparison",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('reference', help='Reference genome FASTA file')
    parser.add_argument('assembly', help='Assembly genome FASTA file')
    parser.add_argument('--kmer-size', type=int, default=201,
                       help='Kmer size (must be odd for canonical representation)')
    parser.add_argument('--bin-size', type=int, default=500000,
                       help='Bin size in base pairs')
    parser.add_argument('--output-prefix', type=str, default='annotation',
                       help='Prefix for output files')
    parser.add_argument('--debug', action='store_true',
                       help='Enable debug logging')
    
    args = parser.parse_args()
    
    if args.debug:
        logger.setLevel(logging.DEBUG)
    
    try:
        # Validate kmer size
        if args.kmer_size < 1:
            raise ValueError("Kmer size must be positive")
        
        # Read input files
        logger.info("Reading reference genome...")
        ref_sequences = read_fasta(args.reference)
        
        logger.info("Reading assembly genome...")
        asm_sequences = read_fasta(args.assembly)
        
        # Process each assembly contig separately
        all_results = []
        
        for asm_id, asm_seq in asm_sequences.items():
            logger.info(f"Processing assembly contig: {asm_id}")
            
            # Extract kmers from this assembly contig
            asm_kmer_list = [extract_kmers_from_sequence(asm_seq, args.kmer_size)]
            asm_id_list = [asm_id]
            
            # Annotate kmers using reference
            kmer_annotation = annotate_kmers(ref_sequences, args.kmer_size, asm_kmer_list)
            
            # Process assembly
            results = process_assembly(
                asm_seq, asm_id, kmer_annotation,
                args.kmer_size, args.bin_size
            )
            all_results.extend(results)
        
        # Write output
        write_output(all_results, args.output_prefix)
        
        logger.info("Annotation completed successfully")
        
    except Exception as e:
        logger.error(f"Error during annotation: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()