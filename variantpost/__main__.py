import argparse
import logging
import math
import multiprocessing as mp
import os
import sys
import tempfile
from collections import namedtuple
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd
import pysam
from scipy.stats import fisher_exact

from .variant import Variant
from .variantalignment import VariantAlignment

# Logging Configuration
logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s - %(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

KEYS_83: Tuple[str, ...] = (
    "1:Del:C:0",
    "1:Del:C:1",
    "1:Del:C:2",
    "1:Del:C:3",
    "1:Del:C:4",
    "1:Del:C:5",
    "1:Del:T:0",
    "1:Del:T:1",
    "1:Del:T:2",
    "1:Del:T:3",
    "1:Del:T:4",
    "1:Del:T:5",
    "1:Ins:C:0",
    "1:Ins:C:1",
    "1:Ins:C:2",
    "1:Ins:C:3",
    "1:Ins:C:4",
    "1:Ins:C:5",
    "1:Ins:T:0",
    "1:Ins:T:1",
    "1:Ins:T:2",
    "1:Ins:T:3",
    "1:Ins:T:4",
    "1:Ins:T:5",
    "2:Del:M:1",
    "2:Del:R:0",
    "2:Del:R:1",
    "2:Del:R:2",
    "2:Del:R:3",
    "2:Del:R:4",
    "2:Del:R:5",
    "2:Ins:R:0",
    "2:Ins:R:1",
    "2:Ins:R:2",
    "2:Ins:R:3",
    "2:Ins:R:4",
    "2:Ins:R:5",
    "3:Del:M:1",
    "3:Del:M:2",
    "3:Del:R:0",
    "3:Del:R:1",
    "3:Del:R:2",
    "3:Del:R:3",
    "3:Del:R:4",
    "3:Del:R:5",
    "3:Ins:R:0",
    "3:Ins:R:1",
    "3:Ins:R:2",
    "3:Ins:R:3",
    "3:Ins:R:4",
    "3:Ins:R:5",
    "4:Del:M:1",
    "4:Del:M:2",
    "4:Del:M:3",
    "4:Del:R:0",
    "4:Del:R:1",
    "4:Del:R:2",
    "4:Del:R:3",
    "4:Del:R:4",
    "4:Del:R:5",
    "4:Ins:R:0",
    "4:Ins:R:1",
    "4:Ins:R:2",
    "4:Ins:R:3",
    "4:Ins:R:4",
    "4:Ins:R:5",
    "5:Del:M:1",
    "5:Del:M:2",
    "5:Del:M:3",
    "5:Del:M:4",
    "5:Del:M:5",
    "5:Del:R:0",
    "5:Del:R:1",
    "5:Del:R:2",
    "5:Del:R:3",
    "5:Del:R:4",
    "5:Del:R:5",
    "5:Ins:R:0",
    "5:Ins:R:1",
    "5:Ins:R:2",
    "5:Ins:R:3",
    "5:Ins:R:4",
    "5:Ins:R:5",
)

KEYS_89: Tuple[str, ...] = (
    "A[Ins(C):R0]A",
    "A[Ins(C):R0]T",
    "Ins(C):R(0,3)",
    "Ins(C):R(4,6)",
    "Ins(C):R(7,9)",
    "A[Ins(T):R(0,4)]A",
    "A[Ins(T):R(0,4)]C",
    "A[Ins(T):R(0,4)]G",
    "C[Ins(T):R(0,4)]A",
    "C[Ins(T):R(0,4)]C",
    "C[Ins(T):R(0,4)]G",
    "G[Ins(T):R(0,4)]A",
    "G[Ins(T):R(0,4)]C",
    "G[Ins(T):R(0,4)]G",
    "A[Ins(T):R(5,7)]A",
    "A[Ins(T):R(5,7)]C",
    "A[Ins(T):R(5,7)]G",
    "C[Ins(T):R(5,7)]A",
    "C[Ins(T):R(5,7)]C",
    "C[Ins(T):R(5,7)]G",
    "G[Ins(T):R(5,7)]A",
    "G[Ins(T):R(5,7)]C",
    "G[Ins(T):R(5,7)]G",
    "A[Ins(T):R(8,9)]A",
    "A[Ins(T):R(8,9)]C",
    "A[Ins(T):R(8,9)]G",
    "C[Ins(T):R(8,9)]A",
    "C[Ins(T):R(8,9)]C",
    "C[Ins(T):R(8,9)]G",
    "G[Ins(T):R(8,9)]A",
    "G[Ins(T):R(8,9)]C",
    "G[Ins(T):R(8,9)]G",
    "Ins(2,4):R0",
    "Ins(5,):R0",
    "Ins(2,4):R1",
    "Ins(5,):R1",
    "Ins(2,):R(2,4)",
    "Ins(2,):R(5,9)",
    "[Del(C):R1]A",
    "[Del(C):R1]T",
    "[Del(C):R2]A",
    "[Del(C):R2]T",
    "[Del(C):R3]A",
    "[Del(C):R3]T",
    "[Del(C):R(4,5)]A",
    "[Del(C):R(4,5)]T",
    "[Del(C):R(1,5)]G",
    "Del(C):R(6,9)",
    "A[Del(T):R(1,4)]A",
    "A[Del(T):R(1,4)]C",
    "A[Del(T):R(1,4)]G",
    "C[Del(T):R(1,4)]A",
    "C[Del(T):R(1,4)]C",
    "C[Del(T):R(1,4)]G",
    "G[Del(T):R(1,4)]A",
    "G[Del(T):R(1,4)]C",
    "G[Del(T):R(1,4)]G",
    "A[Del(T):R(5,7)]A",
    "A[Del(T):R(5,7)]C",
    "A[Del(T):R(5,7)]G",
    "C[Del(T):R(5,7)]A",
    "C[Del(T):R(5,7)]C",
    "C[Del(T):R(5,7)]G",
    "G[Del(T):R(5,7)]A",
    "G[Del(T):R(5,7)]C",
    "G[Del(T):R(5,7)]G",
    "A[Del(T):R(8,9)]A",
    "A[Del(T):R(8,9)]C",
    "A[Del(T):R(8,9)]G",
    "C[Del(T):R(8,9)]A",
    "C[Del(T):R(8,9)]C",
    "C[Del(T):R(8,9)]G",
    "G[Del(T):R(8,9)]A",
    "G[Del(T):R(8,9)]C",
    "G[Del(T):R(8,9)]G",
    "Del(2,4):R1",
    "Del(5,):R1",
    "Del(2,8):U(1,2):R(2,4)",
    "Del(2,):U(1,2):R(5,9)",
    "Del(3,):U(3,):R2",
    "Del(3,):U(3,):R(3,9)",
    "Del(2,5):M1",
    "Del(3,5):M2",
    "Del(4,5):M(3,4)",
    "Del(6,):M1",
    "Del(6,):M2",
    "Del(6,):M3",
    "Del(6,):M(4,)",
    "Complex",
)

ALLOWED_BASES: Set[str] = set("ATCG")

# Global multiprocessing state
_worker_tumor_bam = None
_worker_normal_bam = None
_worker_fasta = None


Result = namedtuple(
    "Result",
    [
        "ComplexPos",
        "ComplexRef",
        "ComplexAlt",
        "TumorSupportingCountFw",
        "TumorSupportingCountRv",
        "TumorNonSupportingCountFw",
        "TumorNonSupportingCountRv",
        "TumorUniqueSupportingCount",
        "TumorUniqueNonSupportingCount",
        "NormalSupportingCountFw",
        "NormalSupportingCountRv",
        "NormalNonSupportingCountFw",
        "NormalNonSupportingCountRv",
        "NormalUniqueSupportingCount",
        "NormalUniqueNonSupportingCount",
        "RefIndelChannelCOSMIC83",
        "RefIndelChannel89",
        "PersonalIndelChannelCOSMIC83",
        "PersonalIndelChannel89",
    ],
)


def validate_file(path: str, label: str) -> None:
    """Validate existence of input files."""
    if not os.path.exists(path):
        logging.error("%s file not found: %s", label, path)
        sys.exit(1)


def get_chrom_order(fasta_path: str, df_chroms: pd.Series) -> List[str]:
    """Retrieve ordered chromosome list prioritized by FASTA header."""
    try:
        with pysam.FastaFile(fasta_path) as fa:
            ref_chroms = [str(c) for c in fa.references]
    except Exception as e:
        logging.warning(
            "Could not read references from FASTA (%s). Falling back to observed chromosomes.",
            e,
        )
        ref_chroms = []

    present_chroms = [str(c) for c in df_chroms.unique().tolist()]
    ordered_chroms = [c for c in ref_chroms if c in present_chroms]
    ordered_chroms.extend([c for c in present_chroms if c not in ordered_chroms])

    return ordered_chroms


def is_valid_indel(ref: str, alt: str) -> bool:
    """Check if the given REF/ALT pair represents a valid indel."""
    len_diff = abs(len(ref) - len(alt))
    if len_diff == 0 or len_diff >= 100:
        return False
    return set(ref).issubset(ALLOWED_BASES) and set(alt).issubset(ALLOWED_BASES)


def flag_soft_overlaps(df: pd.DataFrame, window: int) -> pd.DataFrame:
    if df.empty:
        df["is_soft_overlap"] = False
        return df

    df = df.copy()

    df.sort_values(by=["CHROM", "POS"], inplace=True)
    df.reset_index(drop=True, inplace=True)

    same_chrom_prev = df["CHROM"] == df["CHROM"].shift(1)
    same_chrom_next = df["CHROM"] == df["CHROM"].shift(-1)

    dist_prev = np.where(same_chrom_prev, df["POS"] - df["POS"].shift(1), np.inf)
    dist_next = np.where(same_chrom_next, df["POS"].shift(-1) - df["POS"], np.inf)

    df["is_soft_overlap"] = (dist_prev <= window) | (dist_next <= window)

    return df


def extract_indels_to_dataframe(
    vcf_list: List[str],
    filter_sets: Set[str],
    reference: str,
    window: int,
    exclude_filter_sets: Optional[Set[str]] = None,
) -> Tuple[pd.DataFrame, List[str]]:
    """Extract and aggregate unique indel variants from multiple VCF files."""
    variant_dict: Dict[Tuple[str, int, str, str], Set[str]] = {}
    allow_all = "all" in filter_sets or "." in filter_sets
    exclude_filter_sets = exclude_filter_sets or set()

    for vcf_path in vcf_list:
        validate_file(vcf_path, "VCF")
        vcf_name = os.path.basename(vcf_path)

        with pysam.VariantFile(vcf_path) as vcf:
            for record in vcf:
                rec_filters = set(record.filter.keys()) or {"PASS"}

                if exclude_filter_sets and not rec_filters.isdisjoint(
                    exclude_filter_sets
                ):
                    continue

                if not allow_all and rec_filters.isdisjoint(filter_sets):
                    continue

                ref = record.ref
                for alt in record.alts or []:
                    if alt and is_valid_indel(ref, alt):
                        key = (str(record.chrom), record.pos, ref, alt)
                        variant_dict.setdefault(key, set()).add(vcf_name)

    rows = [
        {
            "CHROM": chrom,
            "POS": pos,
            "REF": ref,
            "ALT": alt,
            "SOURCES": ",".join(sorted(sources)),
            "SUPPORT_COUNT": len(sources),
        }
        for (chrom, pos, ref, alt), sources in variant_dict.items()
    ]

    df = pd.DataFrame(rows)
    chrom_order = []
    if not df.empty:
        df["CHROM"] = df["CHROM"].astype(str)
        chrom_order = get_chrom_order(reference, df["CHROM"])
        chrom_category = pd.CategoricalDtype(categories=chrom_order, ordered=True)
        df["CHROM"] = df["CHROM"].astype(chrom_category)
        df.sort_values(by=["CHROM", "POS"], inplace=True)
        df.reset_index(drop=True, inplace=True)

        if window >= 0:
            df = flag_soft_overlaps(df, window)
            df = df[df["is_soft_overlap"]]
            df.reset_index(drop=True, inplace=True)

    return df, chrom_order


def init_worker(tumor_path: str, normal_path: str, fasta_path: str) -> None:
    """Initialize worker process resources for multiprocessing."""
    global _worker_tumor_bam, _worker_normal_bam, _worker_fasta
    _worker_tumor_bam = pysam.AlignmentFile(tumor_path, "rb")
    _worker_normal_bam = pysam.AlignmentFile(normal_path, "rb")
    _worker_fasta = pysam.FastaFile(fasta_path)


def process_single_row(row: pd.Series) -> Result:
    """Process a single variant record across BAM and FASTA contexts."""
    global _worker_tumor_bam, _worker_normal_bam, _worker_fasta

    if _worker_tumor_bam is None:
        raise RuntimeError("Worker process is uninitialized.")

    v = Variant(str(row["CHROM"]), row["POS"], row["REF"], row["ALT"], _worker_fasta)
    valn = VariantAlignment(v, _worker_tumor_bam, _worker_normal_bam)

    cnt = valn.count_alleles()
    tumor_cnt, normal_cnt = cnt.first, cnt.second
    ref_txn_cosmic, ref_txn_89, personal_txn_cosmic, personal_txn_89 = valn.taxonomize()

    return Result(
        ComplexPos=valn.variant.cpos,
        ComplexRef=valn.variant.cref,
        ComplexAlt=valn.variant.calt,
        TumorSupportingCountFw=tumor_cnt.s_fw,
        TumorSupportingCountRv=tumor_cnt.s_rv,
        TumorNonSupportingCountFw=tumor_cnt.n_fw,
        TumorNonSupportingCountRv=tumor_cnt.n_rv,
        TumorUniqueSupportingCount=tumor_cnt.s,
        TumorUniqueNonSupportingCount=tumor_cnt.n,
        NormalSupportingCountFw=normal_cnt.s_fw,
        NormalSupportingCountRv=normal_cnt.s_rv,
        NormalNonSupportingCountFw=normal_cnt.n_fw,
        NormalNonSupportingCountRv=normal_cnt.n_rv,
        NormalUniqueSupportingCount=normal_cnt.s,
        NormalUniqueNonSupportingCount=normal_cnt.n,
        RefIndelChannelCOSMIC83=ref_txn_cosmic,
        RefIndelChannel89=ref_txn_89,
        PersonalIndelChannelCOSMIC83=personal_txn_cosmic,
        PersonalIndelChannel89=personal_txn_89,
    )


def process_chunk_stream(df_chunk: pd.DataFrame) -> List[Tuple[dict, dict]]:
    """Worker task to process a DataFrame chunk and yield dictionaries."""
    return [
        (row.to_dict(), process_single_row(row)._asdict())
        for _, row in df_chunk.iterrows()
    ]


def filter_row(
    row: pd.Series,
    min_strand_bias_read_count: int,
    strand_bias_thresh: float,
    max_normal_vaf: float,
    min_enrichment_thresh: float,
    pval_thresh: float,
) -> Tuple[float, float, float, str]:
    """Calculate VAF, Strand Bias (SOR), and filter flags for a processed variant."""
    ref_fw, ref_rv = row["TumorNonSupportingCountFw"], row["TumorNonSupportingCountRv"]
    alt_fw, alt_rv = row["TumorSupportingCountFw"], row["TumorSupportingCountRv"]

    t_ref, t_alt = (
        row["TumorUniqueNonSupportingCount"],
        row["TumorUniqueSupportingCount"],
    )
    n_ref, n_alt = (
        row["NormalUniqueNonSupportingCount"],
        row["NormalUniqueSupportingCount"],
    )

    filters = []

    # Strand bias calculation
    sor = 0.0
    if t_alt > min_strand_bias_read_count:
        rf, rr = ref_fw + 1, ref_rv + 1
        af, ar = alt_fw + 1, alt_rv + 1

        r_ref = max(rf / rr, rr / rf)
        r_alt = max(af / ar, ar / af)
        ref_alt_ratio = (af * rr) / (ar * rf)

        if ref_alt_ratio < 1.0:
            ref_alt_ratio = 1.0 / ref_alt_ratio

        try:
            sor = round(math.log(r_ref * r_alt) + math.log(ref_alt_ratio), 4)
        except ValueError:
            sor = 0.0

    if sor >= strand_bias_thresh:
        filters.append("StrandBias")

    t_vaf = t_alt / (t_ref + t_alt) if (t_ref + t_alt) > 0 else 0.0
    n_vaf = n_alt / (n_ref + n_alt) if (n_ref + n_alt) > 0 else 0.0

    if t_vaf == 0:
        filters.append("FailedToDetect")

    if n_vaf > 0:
        if n_vaf >= max_normal_vaf:
            filters.append("PossibleGermline")

        if (t_vaf / n_vaf) < min_enrichment_thresh:
            filters.append("WeakEnrichmentInTumor")

        _, p_val = fisher_exact([[t_alt, t_ref], [n_alt, n_ref]], alternative="greater")
        if p_val > pval_thresh:
            filters.append("WeakStatisticalEvidence")

    filter_str = "PASS" if not filters else ";".join(filters)

    return t_vaf, n_vaf, sor, filter_str


def parse_filter_str(filter_str: Optional[str]) -> Set[str]:
    if not filter_str:
        return set()
    normalized = filter_str.replace(";", ",")
    return {f.strip() for f in normalized.split(",") if f.strip()}


def personalizer(args: argparse.Namespace) -> None:
    """Subcommand handler: Personalize and reanalyze somatic indels."""
    validate_file(args.tumor_bam, "Tumor BAM")
    validate_file(args.normal_bam, "Normal BAM")
    validate_file(args.reference, "Reference FASTA")

    filter_sets = parse_filter_str(args.filters)
    exclude_filter_sets = parse_filter_str(args.exclude_filters)

    logging.info("Extracting indels from VCF files...")
    df, chrom_order = extract_indels_to_dataframe(
        args.vcf, filter_sets, args.reference, args.overlap_window, exclude_filter_sets
    )

    if df.empty:
        logging.warning("No variants to process.")
        return

    num_workers = min(args.processes, mp.cpu_count())
    chunks = [df[i::num_workers] for i in range(num_workers)]

    header_written = False

    logging.info("Processing variants with %d worker(s)...", num_workers)
    with tempfile.TemporaryDirectory(prefix="indelinside_") as tmp_dir:
        tmp_output = os.path.join(tmp_dir, "raw_metrics.tsv")

        with open(tmp_output, "w", encoding="utf-8") as out_f:
            with mp.Pool(
                processes=num_workers,
                initializer=init_worker,
                initargs=(args.tumor_bam, args.normal_bam, args.reference),
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

        logging.info("Applying post-processing filters...")

        dfo = pd.read_csv(tmp_output, sep="\t", dtype={"CHROM": str})

        chrom_order = [str(c) for c in chrom_order]
        missing_chroms = [c for c in dfo["CHROM"].unique() if c not in chrom_order]
        if missing_chroms:
            chrom_order.extend(missing_chroms)

        chrom_category = pd.CategoricalDtype(categories=chrom_order, ordered=True)
        dfo["CHROM"] = dfo["CHROM"].astype(chrom_category)
        dfo.sort_values(by=["CHROM", "POS"], inplace=True)

        res = dfo.apply(
            filter_row,
            min_strand_bias_read_count=args.min_strand_bias_read_count,
            strand_bias_thresh=args.strand_bias_thresh,
            max_normal_vaf=args.max_normal_vaf,
            min_enrichment_thresh=args.tumor_normal_vaf_ratio_thresh,
            pval_thresh=args.fisher_exact_pval_thresh,
            axis=1,
        )

        dfo["TumorVAF"], dfo["NormalVAF"], dfo["StrandBias"], dfo["FILTER"] = zip(*res)
        dfo.to_csv(args.output, sep="\t", index=False)
        logging.info("Successfully wrote results to %s", args.output)


def _generate_matrix(
    df: pd.DataFrame,
    keys_tuple: Tuple[str, ...],
    ref_col: str,
    personal_col: str,
    sample_name: str,
    out_filename: str,
) -> None:
    """Helper method to construct and output mutation spectrum matrices."""
    ref_counts = df[ref_col].value_counts().reindex(keys_tuple, fill_value=0)
    personal_counts = df[personal_col].value_counts().reindex(keys_tuple, fill_value=0)

    matrix_df = pd.DataFrame(
        {
            "MutationType": keys_tuple,
            f"{sample_name}": ref_counts.values,
            f"{sample_name}_personalized": personal_counts.values,
        }
    )

    matrix_df.to_csv(out_filename, sep="\t", index=False)
    logging.info("Generated matrix file: %s", out_filename)


def mut_tbl(args: argparse.Namespace) -> None:
    """Subcommand handler: Build matrix for Indel Signature Analysis."""
    validate_file(args.input, "Input TSV")

    df = pd.read_csv(args.input, sep="\t", dtype={"CHROM": str})

    consensus_level = args.consensus_level
    if consensus_level >= 1:
        df = df[df["SUPPORT_COUNT"] >= consensus_level]
    elif consensus_level == -1:
        df = df[df["SUPPORT_COUNT"] == consensus_level]
    else:
        raise ValueError(f"Invalid consensus level: {consensus_level}")

    # Generate COSMIC 83 Matrix
    _generate_matrix(
        df=df,
        keys_tuple=KEYS_83,
        ref_col="RefIndelChannelCOSMIC83",
        personal_col="PersonalIndelChannelCOSMIC83",
        sample_name=args.sample,
        out_filename=f"{args.sample}_indel_83_matrix.txt",
    )

    # Generate 89 Channel Matrix
    _generate_matrix(
        df=df,
        keys_tuple=KEYS_89,
        ref_col="RefIndelChannel89",
        personal_col="PersonalIndelChannel89",
        sample_name=args.sample,
        out_filename=f"{args.sample}_indel_89_matrix.txt",
    )


def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(prog="indelinside")
    subparsers = parser.add_subparsers(
        dest="subcommand", required=True, help="sub-command help"
    )

    # Subcommand: personalize
    parser_personal = subparsers.add_parser(
        "personalize",
        help="Reanalyze somatic indels on locally personalized reference genome",
    )
    parser_personal.add_argument(
        "-t", "--tumor_bam", required=True, help="Path to tumor BAM file"
    )
    parser_personal.add_argument(
        "-n", "--normal_bam", required=True, help="Path to normal BAM file"
    )
    parser_personal.add_argument(
        "-r", "--reference", required=True, help="Path to reference FASTA file"
    )
    parser_personal.add_argument(
        "-o", "--output", required=True, type=str, help="Path to output TSV file"
    )
    parser_personal.add_argument(
        "-v",
        "--vcf",
        nargs="+",
        required=True,
        help="Path(s) to VCF(s) (space-separated)",
    )
    parser_personal.add_argument(
        "--overlap-window",
        type=int,
        default=-1,
        help=(
            "Window size in bp to detect soft-overlapping indels. "
            "If set (>= 0), filters the output to ONLY include soft-overlapping indels "
            "within POS +/- window. (default: -1, process all indels)"
        ),
    )
    parser_personal.add_argument(
        "--filters",
        type=str,
        default="PASS",
        help="Commna-separated(or semicolon-separated) list of accepted FILTER values (default: PASS). Use 'all' or '.' to allow any.",
    )
    parser_personal.add_argument(
        "--exclude-filters",
        type=str,
        default=None,
        help="Comma-separated (or semicolon-separated) list of FILTER values to exclude (e.g. 'LowQuality,Germline').",
    )
    parser_personal.add_argument(
        "-p",
        "--processes",
        type=int,
        default=4,
        help="Number of parallel worker processes (default: 4)",
    )
    parser_personal.add_argument(
        "--max_normal_vaf", type=float, default=0.1, help="Maximum allowed Normal VAF"
    )
    parser_personal.add_argument(
        "--tumor_normal_vaf_ratio_thresh",
        type=float,
        default=5.0,
        help="Minimum Tumor/Normal VAF enrichment ratio",
    )
    parser_personal.add_argument(
        "--strand_bias_thresh",
        type=float,
        default=7.0,
        help="Strand bias threshold (SOR)",
    )
    parser_personal.add_argument(
        "--min_strand_bias_read_count",
        type=int,
        default=3,
        help="Minimum supporting reads required to calculate SOR",
    )
    parser_personal.add_argument(
        "--fisher_exact_pval_thresh",
        type=float,
        default=0.05,
        help="Fisher's exact test p-value threshold",
    )
    parser_personal.set_defaults(handler=personalizer)

    # Subcommand: matrix
    parser_table = subparsers.add_parser(
        "matrix", help="Prepare mutation matrix for indel signature analysis"
    )
    parser_table.add_argument(
        "-i",
        "--input",
        required=True,
        type=str,
        help="Path to output TSV file from personalize subcommand",
    )
    parser_table.add_argument(
        "-s", "--sample", required=True, type=str, help="Sample identifier string"
    )
    parser_table.add_argument(
        "-c",
        "--consensus_level",
        required=True,
        type=int,
        help="The number of variant callers (consensus level) detected the indel. For example, indels detected by 2 or more callers will be included for analysis if 2 is specified.",
    )
    parser_table.add_argument(
        "--filters",
        type=str,
        default="PASS",
        help="Semicolon-separated list of accepted FILTER values (default: PASS). Use 'all' or '.' to allow any.",
    )
    parser_table.set_defaults(handler=mut_tbl)

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    return parser.parse_args()


def main() -> None:
    args = parse_arguments()
    args.handler(args)


if __name__ == "__main__":
    main()
