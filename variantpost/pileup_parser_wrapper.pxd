from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool as bool_t

cdef void  test_it(
            string &,
            int,
            string &,
            string &,
            int, 
            int,
            string &, 
            vector[string] &,
            vector[bool_t] &, 
            vector[string] &,
            vector[int] &, 
            vector[int] &,
            vector[string] &,
            vector[string] &,           
            vector[vector[int]] &, 
            vector[int] & 
)

