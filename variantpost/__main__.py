#!/usr/bin/env python3
import os
import sys
import pysam
import argparse
import pandas as pd

from .variant import Variant
from .variantalignment import VariantAlignment

ALLOWED_BASES = set("ATCG")
_worker_tumor_bam = None
_worker_normal_bam = None
_worker_fasta = None


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
    parser.add_argument(
        "-v", "--vcf", nargs="+", help="path(s) to VCF(s) (space-separated)"
    )
    parser.add_argument("--vcf_list", help="text file containing paths to VCF")

    parser.add_argument(
        "--filters",
        type=str,
        default="PASS",
        help="Semicolon-separated list of accepted FILTER values (default: PASS). Use 'all' or '.' to allow any.",
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
    vcf_path = None
    if args.vcf:
        vcf_path = args.vcf
    elif args.vcf_list:
        paths = load_paths_from_file(args.vcf_list)
        if not paths:
            sys.stderr.write(f"Error: No paths found in {args.vcf_list}\n")
            sys.exit(1)
        vcf_path = paths
    else:
        sys.stderr.write("Error: Either -v/--vcf or --vcf_list is required.\n\n")
        parser.print_help()
        sys.exit(1)

    validate_file(t_bam, "Tumor BAM")
    validate_file(n_bam, "Normal BAM")

    if isinstance(vcf_path, list):
        for path in vcf_path:
            validate_file(path, "VCF")
    else:
        validate_file(vcf_path, "VCF")

    allowed_filters = {f.strip() for f in args.filters.split(";") if f.strip()}

    return t_bam, n_bam, reference, vcf_path, allowed_filters


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

    allow_all = "all" in filter_set or "." in filter_set

    for vcf_path in vcf_list:
        vcf_name = os.path.basename(vcf_path)

        with pysam.VariantFile(vcf_path) as vcf:
            for record in vcf:
                rec_filters = list(record.filter.keys())
                if not rec_filters:
                    rec_filters = ["PASS"]

                if not allow_all and not filter_set.intersection(rec_filters):
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

    # ----------------------------------------------------
    # TODO: dd dddd
    result1 = "something1"
    result2 = "something2"
    # ----------------------------------------------------

    return result1, result2


def process_chunk(df_chunk):
    results = []
    for _, row in df_chunk.iterrows():
        r1, r2 = process_single_row(row)
        results.append((r1, r2))
    return results


def parallel_process_dataframe(df, tumor_bam, normal_bam, fasta, num_workers=None):
    if df.empty:
        return [], []

    if num_workers is None:
        num_workers = max(1, mp.cpu_count() - 1)

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

    result1_list, result2_list = zip(*all_results)
    return list(result1_list), list(result2_list)


def main():
    tumor_bam, normal_bam, vcf_list, filter_sets = parse_arguments()
    df = extract_indels_to_dataframe(vcf_list, filter_sets)

    if not df.empty:
        r1, r2 = parallel_process_dataframe(
            df=df, tumor_bam=tumor_bam, normal_bam=normal_bam, fasta=reference
        )
        df["result1"] = r1
        df["result2"] = r2
    else:
        df["result1"] = []
        df["result2"] = []


if __name__ == "__main__":
    main()
