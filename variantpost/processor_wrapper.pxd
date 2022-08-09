from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool as bool_t

cdef string  test_it(
                string &,
                string &,
                int,
                string &,
                string &,
                int,
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

