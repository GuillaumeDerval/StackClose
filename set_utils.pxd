# distutils: language = c++

from libcpp cimport bool

cdef int compute_intersection(int n, long* a, long* b, long* dest)
cdef bool is_subset(int n, long* a, long* b)
cdef void add_to_set(int elem_to_add, long* the_set)
cdef bool is_present(int elem, long* the_set)
cdef void copy_set(int n, long* a, long* b)
cdef void copy_union_set(int n, long* a, long* b)
cdef void init_set(int n, long* a)
cpdef compute_supports(file)