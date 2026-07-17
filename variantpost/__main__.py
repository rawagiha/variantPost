#!/usr/bin/env python3
import os
import sys
import pysam
import argparse
import pandas as pd
import multiprocessing as mp
from collections import namedtuple

from .variant import Variant
from .variantalignment import VariantAlignment

ALLOWED_BASES = set("ATCG")
_worker_tumor_bam = None
_worker_normal_bam = None
_worker_fasta = None


Result = namedtuple(
    "Result",
    [
        "TumorSuppotingCountFw",
        "TumorSuppotingCountRv",
        "TumorNonSuppotingCountFw",
        "TumorNonSuppotingCountRv",
        "NormalSuppotingCountFw",
        "NormalSuppotingCountRv",
        "NormalNonSuppotingCountFw",
        "NormalNonSuppotingCountRv",
        "RefIndelChannelCOSMIC83",
        "RefIndelChannel89",
        "PersonalIndelChannelCOSMIC83",
        "PersonalIndelChannel89",
    ],
)


def validate_file(path, label):
    if not os.path.exists(path):
        sys.stderr.write(f"Error: {label} file not found: {path}\n")
        sys.exit(1)


def parse_arguments():
    parser = argparse.ArgumentParser(
        prog="indelinside",
    )

    parser.add_argument("-t", "--tumor_bam", help="path to tumor BAM")
    parser.add_argument("-n", "--normal_bam", help="path to normal BAM")
    parser.add_argument("-r", "--reference", help="path to reference Fasta file")
    parser.add_argument("-o", "--output", type=str, help="path to output TSV file")

    vcf_group = parser.add_mutually_exclusive_group(required=True)
    vcf_group.add_argument(
        "-v", "--vcf", nargs="+", help="path(s) to VCF(s) (space-separated)"
    )
    vcf_group.add_argument("--vcf_list", help="text file containing paths to VCF")

    parser.add_argument(
        "--filters",
        type=str,
        default="PASS",
        help="Semicolon-separated list of accepted FILTER values (default: PASS). Use 'all' or '.' to allow any.",
    )

    parser.add_argument(
        "-p",
        "--processes",
        type=int,
        default=4,
        help="Number of parallel worker processes (default: 4)",
    )

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    t_bam = None
    if args.tumor_bam:
        t_bam = args.tumor_bam
    else:
        sys.stderr.write("Error: -t/--tumor_bam is required.\n\n")
        parser.print_help()
        sys.exit(1)

    n_bam = None
    if args.normal_bam:
        n_bam = args.normal_bam
    else:
        sys.stderr.write("Error: -n/--normal_bam is required.\n\n")
        parser.print_help()
        sys.exit(1)

    reference = None
    if args.reference:
        reference = args.reference
    else:
        sys.stderr.write("Error: -r/--reference is required.\n\n")
        parser.print_help()
        sys.exit(1)

    vcf_path = None
    if args.vcf:
        vcf_path = args.vcf
    else:
        vcf_paths = load_paths_from_file(args.vcf_list)
        if not vcf_paths:
            sys.stderr.write(f"Error: No paths found in {args.vcf_list}\n")
            sys.exit(1)

    validate_file(t_bam, "Tumor BAM")
    validate_file(n_bam, "Normal BAM")

    if isinstance(vcf_path, list):
        for path in vcf_path:
            validate_file(path, "VCF")
    else:
        validate_file(vcf_path, "VCF")

    allowed_filters = {f.strip() for f in args.filters.split(";") if f.strip()}

    return (
        t_bam,
        n_bam,
        reference,
        args.output,
        vcf_path,
        allowed_filters,
        args.processes,
    )


def is_valid_indel(ref, alt):
    if len(ref) == len(alt):
        return False

    if abs(len(ref) - len(alt)) >= 100:
        return False

    if not set(ref).issubset(ALLOWED_BASES) or not set(alt).issubset(ALLOWED_BASES):
        return False

    return True


def extract_indels_to_dataframe(vcf_list, filter_sets):
    # key: (CHROM, POS, REF, ALT) -> value: set of source VCF paths
    variant_dict = {}

    allow_all = "all" in filter_sets or "." in filter_sets

    for vcf_path in vcf_list:
        vcf_name = os.path.basename(vcf_path)

        with pysam.VariantFile(vcf_path) as vcf:
            for record in vcf:
                rec_filters = list(record.filter.keys())
                if not rec_filters:
                    rec_filters = ["PASS"]

                if not allow_all and not filter_sets.intersection(rec_filters):
                    continue

                for alt in record.alts:
                    if alt is None:
                        continue

                    ref = record.ref

                    if is_valid_indel(ref, alt):
                        key = (record.chrom, record.pos, ref, alt)

                        if key not in variant_dict:
                            variant_dict[key] = set()
                        variant_dict[key].add(vcf_name)

    rows = []
    for (chrom, pos, ref, alt), sources in variant_dict.items():
        rows.append(
            {
                "CHROM": chrom,
                "POS": pos,
                "REF": ref,
                "ALT": alt,
                "SOURCES": ",".join(sorted(sources)),
                "SUPPORT_COUNT": len(sources),
            }
        )

    df = pd.DataFrame(rows)

    if not df.empty:
        df = df.sort_values(by=["CHROM", "POS"]).reset_index(drop=True)

    return df


def init_worker(tumor_path, normal_path, fasta_path):
    global _worker_tumor_bam, _worker_normal_bam, _worker_fasta
    _worker_tumor_bam = pysam.AlignmentFile(tumor_path, "rb")
    _worker_normal_bam = pysam.AlignmentFile(normal_path, "rb")
    _worker_fasta = pysam.FastaFile(fasta_path)


def process_single_row(row):
    global _worker_tumor_bam, _worker_normal_bam, _worker_fasta

    if _worker_tumor_bam is None:
        raise RuntimeError("Worker is not initialized with BAM/Fasta files.")

    chrom = row["CHROM"]
    pos = row["POS"]
    ref = row["REF"]
    alt = row["ALT"]

    v = Variant(chrom, pos, ref, alt, _worker_fasta)

    valn = VariantAlignment(v, _worker_tumor_bam, _worker_normal_bam)

    cnt = valn.count_alleles()
    tumor_cnt, normal_cnt = cnt.first, cnt.second

    taxon = valn.taxonomize()
    ref_txn_cosmic, ref_txn_89, personal_txn_cosmic, personal_txn_89 = taxon

    return Result(
        TumorSuppotingCountFw=tumor_cnt.s_fw,
        TumorSuppotingCountRv=tumor_cnt.s_rv,
        TumorNonSuppotingCountFw=tumor_cnt.n_fw,
        TumorNonSuppotingCountRv=tumor_cnt.n_rv,
        NormalSuppotingCountFw=normal_cnt.s_fw,
        NormalSuppotingCountRv=normal_cnt.s_rv,
        NormalNonSuppotingCountFw=normal_cnt.n_fw,
        NormalNonSuppotingCountRv=normal_cnt.n_rv,
        RefIndelChannelCOSMIC83=ref_txn_cosmic,
        RefIndelChannel89=ref_txn_89,
        PersonalIndelChannelCOSMIC83=personal_txn_cosmic,
        PersonalIndelChannel89=personal_txn_89,
    )


def process_chunk(df_chunk):
    results = []
    for _, row in df_chunk.iterrows():
        result_obj = process_single_row(row)
        results.append(result_obj)
    return results


def parallel_process_dataframe(df, tumor_bam, normal_bam, fasta, num_workers=None):
    if df.empty:
        return [], []

    if num_workers is None:
        num_workers = min(8, max(1, mp.cpu_count() - 1))

    chunks = [df[i::num_workers] for i in range(num_workers)]

    with mp.Pool(
        processes=num_workers,
        initializer=init_worker,
        initargs=(tumor_bam, normal_bam, fasta),
    ) as pool:
        chunk_results = pool.map(process_chunk, chunks)

    all_results = [None] * len(df)
    for idx_chunk, chunk in enumerate(chunks):
        for idx_row_in_chunk, orig_idx in enumerate(chunk.index):
            all_results[orig_idx] = chunk_results[idx_chunk][idx_row_in_chunk]

    decomposed_results = zip(*all_results)
    return [list(r) for r in decomposed_results]


def process_chunk_stream(df_chunk):
    chunk_results = []

    for _, row in df_chunk.iterrows():
        result_obj = process_single_row(row)

        chunk_results.append((row.to_dict(), result_obj._asdict()))

    return chunk_results


def main():
    (
        tumor_bam,
        normal_bam,
        reference,
        output,
        vcf_list,
        filter_sets,
        processes,
    ) = parse_arguments()

    df = extract_indels_to_dataframe(vcf_list, filter_sets)
    if df.empty:
        print("No variants to process.", file=sys.stderr)
        return

    num_workers = processes
    chunks = [df[i::num_workers] for i in range(num_workers)]

    out_f = open(output, "w", encoding="utf-8")

    header_written = False

    try:
        with mp.Pool(
            processes=num_workers,
            initializer=init_worker,
            initargs=(tumor_bam, normal_bam, reference),
        ) as pool:
            results_generator = pool.imap_unordered(
                process_chunk_stream, chunks, chunksize=1
            )

            for chunk_list in results_generator:
                for orig_row_dict, metrics_dict in chunk_list:
                    output_line_dict = {**orig_row_dict, **metrics_dict}

                    if not header_written:
                        headers = list(output_line_dict.keys())
                        out_f.write("\t".join(headers) + "\n")
                        header_written = True

                    values = [str(output_line_dict[h]) for h in headers]
                    out_f.write("\t".join(values) + "\n")

    finally:
        if output and out_f is not sys.stdout:
            out_f.close()


if __name__ == "__main__":
    main()
