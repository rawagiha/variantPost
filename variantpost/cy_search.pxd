from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool as bool_t

#from pysam.libcalignmentfile cimport AlignmentFile
cdef object  search_target(
                object, 
                object,
                int, #chrom_len
                bint, #exclude dups
                int, #window
                int, #downsample
                string,
                str,
                int,
                string,
                string,
                int, #map qual
                int, #base qual 
                float, #low qual rate
                int, 
                int, 
                int, 
                int,
                int,
                int, # kmer
                int, 
                int,
)

