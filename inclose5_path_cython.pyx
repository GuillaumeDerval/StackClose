# distutils: language = c++
from libc.math cimport ceil

import numpy as np
cimport numpy as np
cimport cython

from cpython.mem cimport PyMem_Malloc, PyMem_Free
from libcpp cimport bool

cdef struct stack_element:
    long* rows
    long* cols
    long* ignored_cols
    int current_column

cdef inline void init_stack_element(stack_element* elem, int n_row_long, int n_col_long):
    elem.rows = <long*>PyMem_Malloc(sizeof(long) * n_row_long)
    elem.cols = <long*>PyMem_Malloc(sizeof(long) * n_col_long)
    elem.ignored_cols = <long*>PyMem_Malloc(sizeof(long) * n_col_long)

cdef inline void free_stack_element(stack_element* elem):
    PyMem_Free(elem.rows)
    PyMem_Free(elem.cols)
    PyMem_Free(elem.ignored_cols)

cdef inline (bool, bool) compute_intersection(int n, long* a, long* b, long* dest):
    """ Computes dest = a & b """
    cdef bool is_equivalent_to_a = True
    cdef bool is_empty = True

    for i in range(n):
        dest[i] = a[i] & b[i]
        is_equivalent_to_a &= (dest[i] == a[i])
        is_empty &= (dest[i] == 0)

    return is_equivalent_to_a, is_empty

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
cdef inline long[:,:] compute_supports(double[:,:] matrix, int n_rows, int n_cols, int n_rows_long):
    cdef long[:,:] supports = np.zeros((n_cols, n_rows_long), dtype=np.long, order="C")

    for i in range(n_rows):
        for j in range(n_cols):
            if matrix[i,j]:
                add_to_set(i, &supports[j,0])

    return supports
#
# cdef cset_to_set(int n, long* the_set):
#     s = set()
#     for i in range(n*SIZE):
#         if is_present(i, the_set):
#             s.add(i)
#     return s

def inclose5path(matrix, process):
    return __inclose5path(matrix, process)

def inclose5(matrix, process):
    return __inclose5(matrix, process)

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
cdef __inclose5path(double[:,:] matrix, process):
    cdef unsigned long NB_CHECKS = 0
    cdef unsigned long NB_SUB_CHECKS = 0
    cdef unsigned long NB_COL_CHECKS = 0

    cdef int n_rows = matrix.shape[0]
    cdef int n_cols = matrix.shape[1]
    cdef int n_row_long = <int>ceil((<float>n_rows) / (8 * sizeof(long)))
    cdef int n_col_long = <int>ceil((<float>n_cols) / (8 * sizeof(long)))

    cdef stack_element* stack = <stack_element* >PyMem_Malloc(sizeof(stack_element) * matrix.shape[1])
    for i in range(n_cols):
        init_stack_element(&stack[i], n_row_long, n_col_long)

    cdef long* row_buffer = <long*>PyMem_Malloc(sizeof(long) * n_row_long)

    cdef long[:,:] supports = compute_supports(matrix, n_rows, n_cols, n_row_long)

    cdef int stack_length = 0

    # Add first element
    init_set(n_row_long, stack[0].rows)
    init_set(n_col_long, stack[0].cols)
    init_set(n_col_long, stack[0].ignored_cols)
    stack[0].current_column = 0
    for i in range(n_rows):
        add_to_set(i, stack[0].rows)
    stack_length += 1


    cdef long* cur_rows = <long*>PyMem_Malloc(sizeof(long) * n_row_long)
    cdef long* cur_cols = <long*>PyMem_Malloc(sizeof(long) * n_col_long)
    cdef long* cur_ignored_cols = <long*>PyMem_Malloc(sizeof(long) * n_col_long)

    cdef int cur_first_col
    cdef int stack_length_new
    cdef int subset_found_at = -1
    cdef int j
    
    while stack_length != 0:
        # Copy state
        stack_length -= 1
        copy_set(n_row_long, stack[stack_length].rows, cur_rows) # R
        copy_set(n_col_long, stack[stack_length].cols, cur_cols) # C
        copy_set(n_col_long, stack[stack_length].ignored_cols, cur_ignored_cols) # C union P union N
        cur_first_col = stack[stack_length].current_column # y
        stack_length_new = stack_length

        j = cur_first_col
        while j != n_cols:
            if not is_present(j, cur_ignored_cols):
                NB_COL_CHECKS += 1

                is_equivalent, is_empty = compute_intersection(n_row_long, cur_rows, &supports[j,0], row_buffer)

                if is_empty:
                    add_to_set(j, cur_ignored_cols)
                elif is_equivalent:
                    add_to_set(j, cur_cols)
                    add_to_set(j, cur_ignored_cols)
                else:
                    NB_CHECKS += 1
                    idx_found = 0
                    found_one = False
                    subset_found_at = -1

                    for i in range(stack_length_new):
                        # we should prove that
                        # assert not is_present(stack[i].current_column-1, cur_cols)
                        # if this is always true, it proves that we do always less computation (<=) than the original

                        NB_SUB_CHECKS += 1
                        if is_subset(n_row_long, row_buffer, stack[i].rows):
                            #assert not is_present(stack[i].current_column-1, cur_cols)
                            found_one = True
                            subset_found_at = i
                            break

                    if not found_one:
                        # canonic, add it to the stack
                        copy_set(n_row_long, row_buffer, stack[stack_length_new].rows)
                        init_set(n_col_long, stack[stack_length_new].ignored_cols)
                        stack[stack_length_new].current_column = j+1
                        stack_length_new += 1
                    else:
                        for i in range(max(stack_length, subset_found_at+1), stack_length_new):
                            add_to_set(j, stack[i].ignored_cols)
            j += 1

        process()

        for i in range(stack_length, stack_length_new):
            copy_set(n_col_long, cur_cols, stack[i].cols)
            copy_union_set(n_col_long, cur_ignored_cols, stack[i].ignored_cols)
            add_to_set(stack[i].current_column - 1, stack[i].cols)
            #add_to_set(stack[i].current_column - 1, stack[i].ignored_cols) #not needed

        stack_length = stack_length_new

    for i in range(n_cols):
        free_stack_element(&stack[i])
    PyMem_Free(stack)
    PyMem_Free(row_buffer)
    PyMem_Free(cur_rows)
    PyMem_Free(cur_cols)
    PyMem_Free(cur_ignored_cols)

    return NB_CHECKS, NB_SUB_CHECKS, NB_COL_CHECKS














@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
cdef __inclose5(double[:,:] matrix, process):
    cdef unsigned long NB_CHECKS = 0
    cdef unsigned long NB_SUB_CHECKS = 0
    cdef unsigned long NB_COL_CHECKS = 0

    cdef int n_rows = matrix.shape[0]
    cdef int n_cols = matrix.shape[1]
    cdef int n_row_long = <int>ceil((<float>n_rows) / (8 * sizeof(long)))
    cdef int n_col_long = <int>ceil((<float>n_cols) / (8 * sizeof(long)))

    cdef stack_element* stack = <stack_element* >PyMem_Malloc(sizeof(stack_element) * matrix.shape[1])
    for i in range(n_cols):
        init_stack_element(&stack[i], n_row_long, n_col_long)

    cdef long* row_buffer = <long*>PyMem_Malloc(sizeof(long) * n_row_long)

    cdef long[:,:] supports = compute_supports(matrix, n_rows, n_cols, n_row_long)

    cdef int stack_length = 0

    # Add first element
    init_set(n_row_long, stack[0].rows)
    init_set(n_col_long, stack[0].cols)
    init_set(n_col_long, stack[0].ignored_cols)
    stack[0].current_column = 0
    for i in range(n_rows):
        add_to_set(i, stack[0].rows)
    stack_length += 1


    cdef long* cur_rows = <long*>PyMem_Malloc(sizeof(long) * n_row_long)
    cdef long* cur_cols = <long*>PyMem_Malloc(sizeof(long) * n_col_long)
    cdef long* cur_ignored_cols = <long*>PyMem_Malloc(sizeof(long) * n_col_long)

    cdef int cur_first_col
    cdef int stack_length_new
    cdef int subset_found_at = -1
    cdef int j

    while stack_length != 0:
        # Copy state
        stack_length -= 1
        copy_set(n_row_long, stack[stack_length].rows, cur_rows) # R
        copy_set(n_col_long, stack[stack_length].cols, cur_cols) # C
        copy_set(n_col_long, stack[stack_length].ignored_cols, cur_ignored_cols) # C union P union N
        cur_first_col = stack[stack_length].current_column # y
        stack_length_new = stack_length

        j = cur_first_col
        while j != n_cols:
            if not is_present(j, cur_ignored_cols):
                NB_COL_CHECKS += 1

                is_equivalent, is_empty = compute_intersection(n_row_long, cur_rows, &supports[j,0], row_buffer)

                if is_empty:
                    add_to_set(j, cur_ignored_cols)
                elif is_equivalent:
                    add_to_set(j, cur_cols)
                    add_to_set(j, cur_ignored_cols)
                else:
                    NB_CHECKS += 1
                    idx_found = 0
                    found_one = False
                    subset_found_at = -1

                    for i in range(j):
                        if not is_present(i, cur_cols):
                            NB_SUB_CHECKS += 1
                            if is_subset(n_row_long, row_buffer, &supports[i,0]):
                                found_one = True
                                subset_found_at = i
                                break


                    if not found_one:
                        # canonic, add it to the stack
                        copy_set(n_row_long, row_buffer, stack[stack_length_new].rows)
                        init_set(n_col_long, stack[stack_length_new].ignored_cols)
                        stack[stack_length_new].current_column = j+1
                        stack_length_new += 1
                    else:
                        if subset_found_at < cur_first_col:
                            add_to_set(j, cur_ignored_cols)
            j += 1

        process()

        for i in range(stack_length, stack_length_new):
            copy_set(n_col_long, cur_cols, stack[i].cols)
            copy_union_set(n_col_long, cur_ignored_cols, stack[i].ignored_cols)
            add_to_set(stack[i].current_column - 1, stack[i].cols)
            #add_to_set(stack[i].current_column - 1, stack[i].ignored_cols) #not needed

        stack_length = stack_length_new

    for i in range(n_cols):
        free_stack_element(&stack[i])
    PyMem_Free(stack)
    PyMem_Free(row_buffer)
    PyMem_Free(cur_rows)
    PyMem_Free(cur_cols)
    PyMem_Free(cur_ignored_cols)

    return NB_CHECKS, NB_SUB_CHECKS, NB_COL_CHECKS
