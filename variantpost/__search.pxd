from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool as bool_t

cpdef object search_target(
    object bam,
    object second_bam, 
    int chrom_len,
    bint exclude_duplicates,
    int window,
    string reference_file_name,
    str chrom,
    str bam_chrom,
    int pos,
    string ref,
    string alt,
    int mapping_quality_threshold,
    int base_quality_threshold, 
    float low_quality_base_rate_threshold,
    int downsample_threshold,
    int match_score,
    int mismatch_penalty,
    int gap_open_penalty,
    int gap_extension_penalty,
    int kmer_size,
    int local_threshold,
    int retarget_thresh,
    int unspliced_local_reference_start,
    int unspliced_local_reference_start, 
)

