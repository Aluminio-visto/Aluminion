#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Identification, resolution and biological validation of plasmid multimers.

Given an assembly FASTA (and optionally the long-read FASTQ used to build it),
this script detects circular contigs that were assembled as N-mers of the same
underlying replicon, reconstructs a single consensus monomer per multimer, and
optionally validates the call against the raw reads by mapping them back.

Outputs are written to ``--out_path``:

    <basename>_corr.fasta     # input FASTA with multimers collapsed to monomers
    <basename>_report.txt     # statistical report (TODO: re-enable)
    monomers/                 # MAFFT alignments and EMBOSS consensus per monomer
    blastn/                   # per-contig BLAST DBs and self-search outputs
    multimer_mapping/         # minimap2 + samtools artefacts for validation
"""

import argparse
import gzip
import os
import subprocess
import sys

import matplotlib.pyplot as plt
import pandas as pd
from Bio import Align, SeqIO

from _log import get_logger

log = get_logger(__name__)


# =============================================================================
# FASTA / FASTQ helpers
# =============================================================================
def read_fasta(fasta_file, max_len=200000, query_len=1200):
    """Read a FASTA into a ``{contig_id: [sequence, query_seed]}`` mapping.

    Contigs outside the ``(query_len, max_len)`` window are kept with an empty
    seed (they will be passed through unchanged by the downstream pipeline).
    """
    if not os.path.isfile(fasta_file) or os.path.getsize(fasta_file) == 0:
        log.error('Empty FASTA or wrong path: %s', fasta_file)
        sys.exit(1)

    d_fasta = {}
    for seq_record in SeqIO.parse(fasta_file, 'fasta'):
        if query_len < len(seq_record) < max_len:
            query = str(seq_record.seq[:query_len])
            d_fasta[seq_record.id] = [seq_record.seq, query]
        else:
            d_fasta[seq_record.id] = [seq_record.seq, '']
    return d_fasta


def write_fasta(output_file, d_seq):
    """Write a ``{contig_id: sequence}`` mapping as a FASTA file.

    Contig IDs are emitted verbatim — keys must already include the leading
    ``>``. (Kept for backward compatibility with the existing callers.)
    """
    with open(output_file, 'w') as out_handle:
        for header, seq in d_seq.items():
            out_handle.write(f'{header}\n{seq}\n')


# =============================================================================
# BLAST self-search for tandem repeats
# =============================================================================
def create_blastdb(out_path, contig, sequence):
    """Build a single-contig BLAST nucleotide database and return its paths."""
    blast_path = os.path.join(out_path, 'blastn')
    os.makedirs(blast_path, exist_ok=True)
    contig_path = os.path.join(blast_path, f'{contig}.fasta')
    write_fasta(contig_path, {f'>{contig}': sequence})

    log_path = os.path.join(blast_path, f'makeblastdb_{contig}.log')
    cmd = ['makeblastdb', '-in', contig_path, '-dbtype', 'nucl', '-out', contig_path]
    with open(log_path, 'w') as log_handle:
        subprocess.run(cmd, stdout=log_handle, stderr=log_handle, check=True)
    return blast_path, contig_path


def repetition_search(blast_path, contig_path, query_fasta, fasta_basename, contig):
    """Run a BLASTn of the contig seed against itself and return the sorted hits file."""
    blast_out = os.path.join(blast_path, f'{fasta_basename}_{contig}_blast.out')
    blast_log = os.path.join(blast_path, f'blast_{contig}.log')

    blastn_cmd = [
        'blastn', '-query', query_fasta, '-db', contig_path, '-strand', 'plus',
        '-outfmt',
        '6 qseqid sseqid pident qcovhsp length qlen slen qstart qend sstart send sframe evalue bitscore',
        '-perc_identity', '90', '-qcov_hsp_perc', '95',
    ]
    sort_cmd = ['sort', '-n', '-k', '10,11']

    blast = subprocess.run(blastn_cmd, check=True, capture_output=True)
    with open(blast_out, 'w') as std_handle, open(blast_log, 'w') as err_handle:
        subprocess.run(sort_cmd, input=blast.stdout, stdout=std_handle,
                       stderr=err_handle, check=True)
    return blast_out


# =============================================================================
# Monomer extraction and similarity check
# =============================================================================
def sliding_window(pos_list, window=2):
    """Return overlapping windows of size ``window`` over ``pos_list``."""
    return [pos_list[i:i + window] for i in range(len(pos_list) - window + 1)]


def extract_monomers(df_reps, sequence, contig_id):
    """Split a multimer ``sequence`` into candidate monomers using BLAST hit positions.

    Returns ``(monomers_dict, id_of_longest_monomer)``.
    """
    start_positions = list(df_reps[9])
    start_positions.append(len(sequence) + 1)
    positions = sliding_window(start_positions, 2)

    monomers = {}
    largest_len = 0
    largest_id = ''
    for start, end in positions:
        mono_id = f'{contig_id}_{start}_{end - 1}'
        mono_seq = sequence[start - 1:end - 1]
        monomers[mono_id] = mono_seq
        if len(mono_seq) > largest_len:
            largest_len = len(mono_seq)
            largest_id = mono_id

    return monomers, largest_id


def similarity_check(monomers, largest_id):
    """Verify that all candidate monomers are mutually similar to the largest one.

    Returns ``(is_multimer, complete_monomers, partial_monomers)``:

    - ``is_multimer`` is False if any pairwise identity < 0.9.
    - Monomers with identity ≥ 0.9 AND coverage ≥ 0.9 are "complete" (used in
      the multiple alignment); those with identity ≥ 0.9 but coverage < 0.9 are
      "partial" (kept for reporting only).
    """
    candidates = [mid for mid in monomers if mid != largest_id]
    complete = [largest_id]
    partial = []
    is_multimer = True

    for mid in candidates:
        aligner = Align.PairwiseAligner()
        aligner.mode = 'global'
        aln = aligner.align(monomers[largest_id], monomers[mid])[0]

        identity = aln.score / len(monomers[largest_id])
        aligned_query = aln.aligned[0]
        aligned_length = sum(end - start for start, end in aligned_query)
        qcov = aligned_length / len(monomers[largest_id])

        if identity < 0.9:
            is_multimer = False
        elif qcov >= 0.9:
            complete.append(mid)
        else:
            partial.append(mid)

    return is_multimer, complete, partial


# =============================================================================
# Read-length validation
# =============================================================================
def extract_read_lengths_and_filter_reads(fastq_gz_file, monomer_length):
    """Return all read lengths plus the subset of reads longer than 1.1 * monomer length."""
    read_lengths = []
    filtered_reads = []
    threshold = monomer_length * 1.1
    with gzip.open(fastq_gz_file, 'rt') as handle:
        for record in SeqIO.parse(handle, 'fastq'):
            read_lengths.append(len(record.seq))
            if len(record.seq) > threshold:
                filtered_reads.append(record)
    return read_lengths, filtered_reads


def plot_read_lengths(read_lengths, output_file, monomer_length, max_length):
    """Histogram of read lengths with annotated monomer and multimer markers."""
    plt.figure(figsize=(10, 6))
    plt.hist(read_lengths, bins=100, alpha=0.75, edgecolor='black')
    plt.axvline(x=monomer_length, color='red', linestyle='--', linewidth=2)
    plt.text(monomer_length, max(plt.ylim()) * 0.9,
             f'X={monomer_length}', color='red', ha='right')
    plt.axvline(x=max_length, color='blue', linestyle='--', linewidth=2)
    plt.text(max_length, max(plt.ylim()) * 0.8,
             f'X={max_length}', color='blue', ha='right')
    plt.yscale('log')
    plt.title('Distribution of Read Lengths')
    plt.xlabel('Read Length')
    plt.ylabel('Frequency')
    plt.grid(True)
    plt.savefig(output_file)
    plt.close()


def save_filtered_reads(filtered_reads, output_file):
    """Write the subset of long reads (longer than the monomer) back to a gzipped FASTQ."""
    with gzip.open(output_file, 'wt') as handle:
        SeqIO.write(filtered_reads, handle, 'fastq')


# =============================================================================
# CLI
# =============================================================================
def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Identification, resolution and biological confirmation of plasmid multimers.',
    )
    parser.add_argument('--fasta_file', type=str, required=True,
                        help='Path to the plasmid FASTA file.')
    parser.add_argument('--fastq_file', type=str, default=None,
                        help='Path to the long-read FASTQ.gz used to build the assembly.')
    parser.add_argument('--out_path', type=str, default=None,
                        help='Output path (default: <fasta_file>_monomer).')
    parser.add_argument('--threads', type=int, default=6,
                        help='Number of threads (default: 6).')
    parser.add_argument('--max_len', type=int, default=200000,
                        help='Maximum plasmid length to inspect (default: 200000).')

    args = parser.parse_args()
    if args.out_path is None:
        args.out_path = os.path.splitext(args.fasta_file)[0] + '_monomer'
    return args


# =============================================================================
# Main
# =============================================================================
def main():
    args = parse_arguments()

    fasta_file = args.fasta_file
    fasta_basename = os.path.splitext(os.path.basename(fasta_file))[0]
    out_path = args.out_path
    os.makedirs(out_path, exist_ok=True)
    threads = args.threads

    # Read input FASTA and select contigs in the inspection range.
    d_fasta = read_fasta(fasta_file, max_len=args.max_len)

    # -------------------------------------------------------------------------
    # Step 1 — Tandem repeat identification
    # -------------------------------------------------------------------------
    log.info('Step 1: Repetition identification')
    d_contig_corr = {}
    multimers = []
    monomer_len = 0  # Updated when a multimer is resolved (used in step 3).

    for contig in d_fasta:
        log.info('Processing contig: %s', contig)
        sequence = str(d_fasta[contig][0])
        query = d_fasta[contig][1]

        if query == '':
            log.info('  Skipping (contig outside size window).')
            d_contig_corr[f'>{contig}'] = sequence
            continue

        # Build a per-contig BLAST DB and search the seed against the contig.
        blast_path, contig_path = create_blastdb(out_path, contig, sequence)
        query_fasta = os.path.join(blast_path, f'query_{contig}.fasta')
        write_fasta(query_fasta, {'>query_sequence': query})
        blast_out = repetition_search(blast_path, contig_path, query_fasta,
                                      fasta_basename, contig)

        # Extract candidate monomers from the BLAST hit positions.
        df_reps = pd.read_table(blast_out, header=None)
        monomers, largest_id = extract_monomers(df_reps, sequence, contig)
        monomer_len = len(monomers[largest_id])

        if len(monomers) == 1:
            log.info('  Skipping (not a multimer).')
            d_contig_corr[f'>{contig}'] = sequence
            continue

        # Verify all monomers are mutually similar to the longest one.
        is_multimer, l_complete, l_partial = similarity_check(monomers, largest_id)
        if not is_multimer:
            log.info('  Skipping (monomers diverge — not a true multimer).')
            d_contig_corr[f'>{contig}'] = sequence
            continue

        # ---------------------------------------------------------------------
        # Step 2 — Multimer resolution: MAFFT alignment + EMBOSS consensus
        # ---------------------------------------------------------------------
        multimers.append(contig)
        log.info('  Multimer detected — building consensus monomer.')

        monomer_path = os.path.join(out_path, 'monomers')
        os.makedirs(monomer_path, exist_ok=True)
        monomer_file = os.path.join(monomer_path, 'complete_monomers.fasta')
        with open(monomer_file, 'w') as out_handle:
            for mid in l_complete:
                out_handle.write(f'>{mid}\n{monomers[mid]}\n')

        aln_file = os.path.join(monomer_path, 'complete_monomers.aln')
        aln_err = os.path.join(monomer_path, 'mafft.log')
        mafft_cmd = ['mafft', '--adjustdirectionaccurately',
                     '--thread', str(threads), monomer_file]
        with open(aln_file, 'w') as std_handle, open(aln_err, 'w') as err_handle:
            subprocess.run(mafft_cmd, stdout=std_handle, stderr=err_handle, check=True)

        cons_file = os.path.join(out_path, f'{contig}_consensus.fasta')
        cons_cmd = ['em_cons', '-sequence', aln_file, '-outseq', cons_file,
                    '-name', contig]
        subprocess.run(cons_cmd, check=True)
        log.info('  Multimer resolved into a single monomer consensus.')

        d_cons = read_fasta(cons_file, query_len=0)
        d_contig_corr[f'>{contig}'] = str(d_cons[contig][0])

    output_file = os.path.join(out_path, f'{fasta_basename}_corr.fasta')
    write_fasta(output_file, d_contig_corr)
    log.info('Corrected assembly written to: %s', output_file)

    # -------------------------------------------------------------------------
    # Step 3 — Optional biological validation against the raw reads
    # -------------------------------------------------------------------------
    fastq_file = args.fastq_file
    if not multimers or fastq_file is None:
        log.info('Skipping multimer validation (no multimers detected or no FASTQ provided).')
        return

    log.info('Step 3: Biological multimer validation')
    log.info('Mapping reads to assembly')

    map_path = os.path.join(out_path, 'multimer_mapping')
    os.makedirs(map_path, exist_ok=True)

    sam_file = os.path.join(map_path, f'{fasta_basename}_whole_genome.sam')
    sam_err = os.path.join(map_path, 'minimap.log')
    minimap_cmd = ['minimap2', '--secondary=no', '-t', str(threads),
                   '-ax', 'lr:hq', '-o', sam_file, fasta_file, fastq_file]
    with open(sam_err, 'w') as err_handle:
        subprocess.run(minimap_cmd, stderr=err_handle, check=True)

    # SAM → BAM filtering unmapped reads.
    bam_file = os.path.join(map_path, f'{fasta_basename}_whole_genome.bam')
    smt_err = os.path.join(map_path, 'samtools.log')
    samtools_filter_cmd = ['samtools', 'view', '-h', '-@', str(threads),
                           '-F', '4', '-bS', sam_file]
    with open(bam_file, 'w') as std_handle, open(smt_err, 'w') as err_handle:
        subprocess.run(samtools_filter_cmd, stdout=std_handle,
                       stderr=err_handle, check=True)

    # Sort and index.
    bam_sorted_file = os.path.join(map_path, f'{fasta_basename}_whole_genome.sorted.bam')
    samtools_sort_cmd = ['samtools', 'sort', '--threads', str(threads),
                         bam_file, '-o', bam_sorted_file]
    with open(smt_err, 'a') as err_handle:
        subprocess.run(samtools_sort_cmd, stderr=err_handle, check=True)
    with open(smt_err, 'a') as err_handle:
        subprocess.run(['samtools', 'index', '-@', str(threads), bam_sorted_file],
                       stderr=err_handle, check=True)

    # Per-multimer: extract reads mapping the multimer, validate by read-length plot.
    for multimer in multimers:
        log.info('Filtering reads mapping plasmid: %s', multimer)

        sam_concat_file = os.path.join(map_path, f'{fasta_basename}_{multimer}.sam')
        with open(sam_concat_file, 'w') as std_handle, open(smt_err, 'a') as err_handle:
            subprocess.run(['samtools', 'view', '-h', bam_sorted_file, multimer],
                           stdout=std_handle, stderr=err_handle, check=True)

        bam_concat_file = os.path.join(map_path, f'{fasta_basename}_{multimer}.bam')
        with open(bam_concat_file, 'w') as std_handle, open(smt_err, 'a') as err_handle:
            subprocess.run(['samtools', 'view', '-h', '-@', str(threads),
                            '-F', '4', '-bS', sam_concat_file],
                           stdout=std_handle, stderr=err_handle, check=True)

        with open(smt_err, 'a') as err_handle:
            subprocess.run(['samtools', 'index', '-@', str(threads), bam_concat_file],
                           stderr=err_handle, check=True)

        fastq_sorted_file = os.path.join(
            out_path, f'{fasta_basename}_{multimer}.sorted.fastq.gz')
        fastq = subprocess.run(['samtools', 'fastq', bam_concat_file],
                               check=True, capture_output=True)
        with open(fastq_sorted_file, 'w') as std_handle, open(smt_err, 'a') as err_handle:
            subprocess.run(['gzip'], input=fastq.stdout, stdout=std_handle,
                           stderr=err_handle, check=True)

        log.info('Plotting read lengths')
        read_lengths, _filtered_reads = extract_read_lengths_and_filter_reads(
            fastq_sorted_file, monomer_len)
        output_plot = os.path.join(out_path, f'{fasta_basename}_{multimer}_read_plot.png')
        multimer_len = len(d_fasta[multimer][0])
        monomer_len = len(d_contig_corr[f'>{multimer}'])
        plot_read_lengths(read_lengths, output_plot, monomer_len, multimer_len)

        log.info('Plot saved as: %s', output_plot)

    log.info('End of pipeline.')


if __name__ == '__main__':
    main()
