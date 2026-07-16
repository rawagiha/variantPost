#!/usr/bin/env python3
import os
import argparse
import sys

from .variant import Variant
from .variantalignment import VariantAlignment 

def validate_file(path, label):
    if not os.path.exists(path):
        sys.stderr.write(f"Error: {label} file not found: {path}\n")
        sys.exit(1) 

def parse_arguments():
    parser = argparse.ArgumentParser(
        prog="indelinside",
        description="Run variant alignment and somatic analysis."
    )

    parser.add_argument("-t", "--tumor_bam", help="Path to Tumor BAM file")
    parser.add_argument("-n", "--normal_bam", help="Path to Normal BAM file")
    parser.add_argument("-v", "--vcf", nargs="+", help="Path(s) to VCF file(s) (space-separated)")
    parser.add_argument("--vcf_list", help="Text file containing path to VCF")

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

    return t_bam, n_bam, vcf_path

def main():
    tumor_bam, normal_bam, vcf_list = parse_arguments()
    print("indelinside test 2")


if __name__ == "__main__":
    main()
