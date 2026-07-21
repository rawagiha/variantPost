#!/usr/bin/env python3
import os
import sys
import math
import pysam
import argparse
import tempfile
import pandas as pd
import multiprocessing as mp
from scipy.stats import fisher_exact
from collections import namedtuple, OrderedDict


from .variant import Variant
from .variantalignment import VariantAlignment

KEYS_83 = (
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

KEYS_89 = (
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

ALLOWED_BASES = set("ATCG")
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


def validate_file(path, label):
    if not os.path.exists(path):
        sys.stderr.write(f"Error: {label} file not found: {path}\n")
        sys.exit(1)


def get_chrom_order(fasta_path, df_chroms):
    try:
        with pysam.FastaFile(fasta_path) as fa:
            ref_chroms = list(fa.references)
    except Exception:
        ref_chroms = []

    present_chroms = df_chroms.unique().tolist()

    ordered_chroms = [c for c in ref_chroms if c in present_chroms]

    for c in present_chroms:
        if c not in ordered_chroms:
            ordered_chroms.append(c)

    return ordered_chroms


def is_valid_indel(ref, alt):
    if len(ref) == len(alt):
        return False

    if abs(len(ref) - len(alt)) >= 100:
        return False

    if not set(ref).issubset(ALLOWED_BASES) or not set(alt).issubset(ALLOWED_BASES):
        return False

    return True


def extract_indels_to_dataframe(vcf_list, filter_sets, reference):
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
    chrom_order = []
    if not df.empty:
        chrom_order = get_chrom_order(reference, df["CHROM"])
        chrom_category = pd.CategoricalDtype(categories=chrom_order, ordered=True)
        df["CHROM"] = df["CHROM"].astype(chrom_category)
        df.sort_values(by=["CHROM", "POS"], inplace=True)
        df.reset_index(drop=True, inplace=True)

    return df, chrom_order


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


def process_chunk(df_chunk):
    results = []
    for _, row in df_chunk.iterrows():
        result_obj = process_single_row(row)
        results.append(result_obj)
    return results


def parallel_process_dataframe(df, tumor_bam, normal_bam, fasta, num_workers=None):
    if df.empty:
        return [], []

    # df is re-indexed with drop

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


def get_chrom_order(fasta_path, df_chroms):
    try:
        with pysam.FastaFile(fasta_path) as fa:
            ref_chroms = list(fa.references)
    except Exception:
        ref_chroms = []

    present_chroms = df_chroms.unique().tolist()

    ordered_chroms = [c for c in ref_chroms if c in present_chroms]

    for c in present_chroms:
        if c not in ordered_chroms:
            ordered_chroms.append(c)

    return ordered_chroms


def filter_row(
    row,
    min_strand_bias_read_count,
    strand_bias_thresh,
    max_normal_vaf,
    min_enrichment_thresh,
    pval_thresh,
):
    ref_fw = row["TumorNonSupportingCountFw"]
    ref_rv = row["TumorNonSupportingCountRv"]
    alt_fw = row["TumorSupportingCountFw"]
    alt_rv = row["TumorSupportingCountRv"]

    t_ref = row["TumorUniqueNonSupportingCount"]
    t_alt = row["TumorUniqueSupportingCount"]
    n_ref = row["NormalUniqueNonSupportingCount"]
    n_alt = row["NormalUniqueSupportingCount"]

    filters = []
    # Strand bias odds ratio filter
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
            sor = math.log(r_ref * r_alt) + math.log(ref_alt_ratio)
        except ValueError:
            sor = 0.0

        sor = round(sor, 4)

    if sor >= strand_bias_thresh:
        filters.append("StrandBias")

    t_vaf = t_alt / (t_ref + t_alt) if (t_ref + t_alt) > 0 else 0
    n_vaf = n_alt / (n_ref + n_alt) if (n_ref + n_alt) > 0 else 0

    # Not detected
    if t_vaf == 0:
        filters.append("FailedToDetect")

    # Normal VAF > 0
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


def personalizer(args):
    tumor_bam = args.tumor_bam
    normal_bam = args.normal_bam
    reference = args.reference
    output = args.output
    vcf_list = args.vcf
    processes = args.processes

    filter_sets = {f.strip() for f in args.filters.split(";") if f.strip()}

    df, chrom_order = extract_indels_to_dataframe(vcf_list, filter_sets, reference)
    if df.empty:
        print("No variants to process.", file=sys.stderr)
        return

    num_workers = processes
    chunks = [df[i::num_workers] for i in range(num_workers)]

    header_written = False

    with tempfile.TemporaryDirectory(prefix="indelinside_") as tmp_dir:
        tmp_output = os.path.join(tmp_dir, "raw_metrics.tsv")
        out_f = open(tmp_output, "w", encoding="utf-8")

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
            out_f.close()

        dfo = pd.read_csv(tmp_output, sep="\t")

        chrom_category = pd.CategoricalDtype(categories=chrom_order, ordered=True)
        dfo["CHROM"] = dfo["CHROM"].astype(chrom_category)
        dfo.sort_values(by=["CHROM", "POS"], inplace=True)

        (
            dfo["TumorVAF"],
            dfo["NormalVAF"],
            dfo["StrandBias"],
            dfo["FILTER"],
        ) = zip(
            *dfo.apply(
                filter_row,
                min_strand_bias_read_count=args.min_strand_bias_read_count,
                strand_bias_thresh=args.strand_bias_thresh,
                max_normal_vaf=args.max_normal_vaf,
                min_enrichment_thresh=args.tumor_normal_vaf_ratio_thresh,
                pval_thresh=args.fisher_exact_pval_thresh,
                axis=1,
            )
        )

        dfo.to_csv(output, sep="\t", index=False)


def mut_tbl(args):
    df = pd.read_csv(args.input, sep="\t")

    consensus_level = args.consensus_level
    if consensus_level >= 1:
        df = df[df["SUPPORT_COUNT"] >= consensus_level]
    elif consensus_level == -1:
        df = df[df["SUPPORT_COUNT"] == consensus_level]
    else:
        raise ValueError(f"Invalide consensus level: {consensus_level}")

    df.reset_index(drop=True, inplace=True)

    ref_83 = OrderedDict((key, 0) for key in KEYS_83)
    ref_83_cnt = df["RefIndelChannelCOSMIC83"].value_counts()
    for key, cnt in ref_83_cnt.items():
        if key in KEYS_83:
            ref_83[key] = cnt

    personal_83 = OrderedDict((key, 0) for key in KEYS_83)
    personal_83_cnt = df["PersonalIndelChannelCOSMIC83"].value_counts()
    for key, cnt in personal_83_cnt.items():
        if key in KEYS_83:
            personal_83[key] = cnt

    data_83 = []
    for key in KEYS_83:
        d = {
            "MutationType": key,
            f"{args.sample}": ref_83[key],
            f"{args.sample}_personalized": personal_83[key],
        }
        data_83.append(d)

    df_83 = pd.DataFrame(data_83)
    df_83.to_csv(f"{args.sample}_indel_83_matrix.txt", sep="\t", index=False)

    ref_89 = OrderedDict((key, 0) for key in KEYS_89)
    ref_89_cnt = df["RefIndelChannel89"].value_counts()
    for key, cnt in ref_89_cnt.items():
        if key in KEYS_89:
            ref_89[key] = cnt

    personal_89 = OrderedDict((key, 0) for key in KEYS_89)
    personal_89_cnt = df["PersonalIndelChannel89"].value_counts()
    for key, cnt in personal_89_cnt.items():
        if key in KEYS_89:
            personal_89[key] = cnt

    data_89 = []
    for key in KEYS_89:
        d = {
            "MutationType": key,
            f"{args.sample}": ref_89[key],
            f"{args.sample}_personalized": personal_89[key],
        }
        data_89.append(d)

    df_89 = pd.DataFrame(data_89)
    df_89.to_csv(f"{args.sample}_indel_89_matrix.txt", sep="\t", index=False)


def parse_arguments():
    parser = argparse.ArgumentParser(
        prog="indelinside",
    )

    subparsers = parser.add_subparsers(
        dest="subcommand", required=True, help="sub-command help"
    )

    parser_personal = subparsers.add_parser(
        "personalize",
        help="reanalyze somatic indels on locally personalized reference genome",
    )

    parser_personal.add_argument("-t", "--tumor_bam", help="path to tumor BAM")
    parser_personal.add_argument("-n", "--normal_bam", help="path to normal BAM")
    parser_personal.add_argument(
        "-r", "--reference", help="path to reference Fasta file"
    )
    parser_personal.add_argument(
        "-o", "--output", type=str, help="path to output TSV file"
    )

    # vcf_group = parser_personal.add_mutually_exclusive_group(required=True)
    parser_personal.add_argument(
        "-v", "--vcf", nargs="+", help="path(s) to VCF(s) (space-separated)"
    )
    # vcf_group.add_argument("--vcf_list", help="text file containing paths to VCF")

    parser_personal.add_argument(
        "--filters",
        type=str,
        default="PASS",
        help="Semicolon-separated list of accepted FILTER values (default: PASS). Use 'all' or '.' to allow any.",
    )

    parser_personal.add_argument(
        "-p",
        "--processes",
        type=int,
        default=4,
        help="Number of parallel worker processes (default: 4)",
    )

    parser_personal.add_argument("--max_normal_vaf", type=float, default=0.1, help="aa")
    parser_personal.add_argument(
        "--tumor_normal_vaf_ratio_thresh", type=float, default=5.0, help="aa"
    )
    parser_personal.add_argument(
        "--strand_bias_thresh", type=float, default=7.0, help="applied to count > 1"
    )
    parser_personal.add_argument(
        "--min_strand_bias_read_count", type=int, default=3, help="aa"
    )
    parser_personal.add_argument(
        "--fisher_exact_pval_thresh", type=float, default=0.05, help="aa"
    )
    parser_personal.set_defaults(handler=personalizer)

    parser_table = subparsers.add_parser(
        "matrix", help="prepare mutation matrix for indel signature analysis"
    )
    parser_table.add_argument(
        "-i",
        "--input",
        type=str,
        help="path to the output TSV file from personalize subcommand",
    )
    parser_table.add_argument("-s", "--sample", type=str, help="aa")
    parser_table.add_argument("-c", "--consensus_level", type=int, help="aa")
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

    if len(sys.argv) == 2 and sys.argv[1] in ["personalize", "matrix"]:
        if sys.argv[1] == "personalize":
            parser_personal.print_help()
        elif sys.argv[1] == "matrix":
            parser_table.print_help()
        sys.exit(0)

    return parser.parse_args()


def main():
    args = parse_arguments()

    args.handler(args)


if __name__ == "__main__":
    main()
