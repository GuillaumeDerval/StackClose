# distutils: language = c++
from libc.math cimport ceil

import numpy as np
cimport numpy as np
cimport cython

from libcpp cimport bool

from libc.stdio cimport FILE, fopen, fclose, getline
from libc.stdlib cimport free, strtol
from cython.operator cimport dereference as deref, preincrement as inc

cdef extern int __builtin_popcountl(unsigned long x)
cdef extern from "<vector>" namespace "std":
    cdef cppclass vector[T]:
        cppclass iterator:
            T operator*()
            iterator operator++()
            bint operator==(iterator)
            bint operator!=(iterator)
        vector()
        void push_back(T&)
        T& operator[](int)
        T& at(int)
        iterator begin()
        iterator end()

cdef inline int compute_intersection(int n, long* a, long* b, long* dest):
    """ Computes dest = a & b """
    cdef int length = 0

    for i in range(n):
        dest[i] = a[i] & b[i]
        length += __builtin_popcountl(dest[i])

    return length

cdef inline bool is_subset(int n, long* a, long* b):
    """ Returns (a subset of b) """
    for i in range(n):
        if (a[i] & b[i]) != a[i]:
            return False
    return True

cdef int SIZE = sizeof(long)*8

@cython.cdivision(True)
cdef inline void add_to_set(int elem_to_add, long* the_set):
    the_set[elem_to_add // SIZE] |= 1L << (elem_to_add % SIZE)

@cython.cdivision(True)
cdef inline bool is_present(int elem, long* the_set):
    return (the_set[elem // SIZE] & (1L << (elem % SIZE))) != 0L

cdef inline void copy_set(int n, long* a, long* b):
    """ Copy a to b"""
    for i in range(n):
        b[i] = a[i]

cdef inline void copy_union_set(int n, long* a, long* b):
    """ Copy (a|b) to b"""
    for i in range(n):
        b[i] |= a[i]

cdef inline void init_set(int n, long* a):
    """ Init a set: make it empty"""
    for i in range(n):
        a[i] = 0L

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.cdivision(True)
cpdef inline compute_supports(file):
    cdef FILE * fp;
    cdef char * line = NULL
    cdef char * subline = NULL
    cdef char * sublineOut = NULL
    cdef size_t len = 0
    cdef ssize_t read;

    fp = fopen(file.encode('utf-8'), "r")
    if fp == NULL:
        raise Exception("No file!")

    cdef int idx = 0
    cdef int lastRead = 0
    cdef vector[int] seen
    cdef int maxCol = -1

    while True:
        read = getline(&line, &len, fp)
        if read == -1:
            break

        subline = line

        while True:
            lastRead = strtol(subline, &sublineOut, 10)
            if subline == sublineOut:
                break
            subline = sublineOut
            if lastRead < 0:
                raise Exception("Invalid item id")
            maxCol = max(maxCol, lastRead)
            seen.push_back(lastRead)
        seen.push_back(-1)

        idx += 1

    fclose(fp)
    if line:
        free(line)

    cdef int n_rows = idx
    cdef int n_cols = maxCol + 1

    cdef int n_row_long = <int>ceil((<float>n_rows) / (8 * sizeof(long)))
    cdef int n_col_long = <int>ceil((<float>n_cols) / (8 * sizeof(long)))

    cdef long[:,:] supports = np.zeros((n_cols, n_row_long), dtype=np.long, order="C")

    cdef vector[int].iterator itr = seen.begin()
    idx = 0
    while itr != seen.end():
        lastRead = deref(itr)
        inc(itr)

        if lastRead == -1:
            idx += 1
        else:
            add_to_set(idx, &supports[lastRead,0])

    return supports, n_rows, n_cols, n_row_long, n_col_long

cdef inline cset_to_set(int n, long* the_set):
    s = set()
    for i in range(n*SIZE):
        if is_present(i, the_set):
            s.add(i)
    return s