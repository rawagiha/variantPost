from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool as bool_t

cdef object  search_target(
                string &,
                string &,
                int,
                string &,
                string &,
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
                #string &, 
                vector[string] &,
                vector[bool_t] &, 
                vector[string] &,
                vector[int] &, 
                vector[int] &,
                vector[string] &,
                #vector[string] &,           
                vector[vector[int]] &, 
                vector[int] &,
                vector[bool_t] &
)

