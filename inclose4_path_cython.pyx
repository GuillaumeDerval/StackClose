# distutils: language = c++
from collections import deque
from libc.math cimport ceil

import numpy as np
cimport numpy as np

from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free
from libcpp cimport bool

cdef struct stack_element:
    long* rows
    long* cols
    long* ignored_cols
    int current_column

cdef void init_stack_element(stack_element* elem, int n_row_long, int n_col_long):
    elem.rows = <long*>PyMem_Malloc(sizeof(long) * n_row_long)
    elem.cols = <long*>PyMem_Malloc(sizeof(long) * n_col_long)
    elem.ignored_cols = <long*>PyMem_Malloc(sizeof(long) * n_col_long)

cdef void free_stack_element(stack_element* elem):
    PyMem_Free(elem.rows)
    PyMem_Free(elem.cols)
    PyMem_Free(elem.ignored_cols)

cdef (bool, bool) compute_intersection(int n, long* a, long* b, long* dest):
    """ Computes dest = a & b """
    cdef bool is_equivalent_to_a = True
    cdef bool is_empty = True

    for i in range(n):
        dest[i] = a[i] & b[i]
        is_equivalent_to_a &= (dest[i] == a[i])
        is_empty &= (dest[i] == 0)

    return is_equivalent_to_a, is_empty

cdef bool is_subset(int n, long* a, long* b):
    """ Returns (a subset of b) """
    for i in range(n):
        if (a[i] & b[i]) != a[i]:
            return False
    return True

cdef int SIZE = sizeof(long)*8

cdef void add_to_set(int elem_to_add, long* the_set):
    the_set[elem_to_add // SIZE] |= 1L << (elem_to_add % SIZE)

cdef bool is_present(int elem, long* the_set):
    return (the_set[elem // SIZE] & (1L << (elem % SIZE))) != 0L

cdef void copy_set(int n, long* a, long* b):
    """ Copy a to b"""
    for i in range(n):
        b[i] = a[i]

cdef void init_set(int n, long* a):
    """ Init a set: make it empty"""
    for i in range(n):
        a[i] = 0

cdef long[:,:] compute_supports(double[:,:] matrix, int n_rows, int n_cols, int n_rows_long):
    cdef long[:,:] supports = np.zeros((n_cols, n_rows_long), dtype=np.long, order="C")

    for i in range(n_rows):
        for j in range(n_cols):
            if matrix[i,j]:
                add_to_set(i, &supports[j,0])

    return supports

cdef cset_to_set(int n, long* the_set):
    s = set()
    for i in range(n*SIZE):
        if is_present(i, the_set):
            s.add(i)
    return s

def inclose5path(matrix, process):
    return __inclose5path(matrix, process)

cdef __inclose5path(double[:,:] matrix, process):
    cdef int NB_CHECKS = 0
    cdef int NB_SUB_CHECKS = 0
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
    
    while stack_length != 0:
        # Copy state
        stack_length -= 1
        copy_set(n_row_long, stack[stack_length].rows, cur_rows) # R
        copy_set(n_col_long, stack[stack_length].cols, cur_cols) # C
        copy_set(n_col_long, stack[stack_length].ignored_cols, cur_ignored_cols) # C union P union N
        cur_first_col = stack[stack_length].current_column # y
        stack_length_new = stack_length

        #print(len(cset_to_set(n_row_long, cur_rows)), cur_first_col)
        #print("{:09b} {:09b} {:09b} {}".format(cur_rows[0], cur_cols[0], cur_ignored_cols[0], cur_first_col))

        for j in range(cur_first_col, n_cols):
            if not is_present(j, cur_cols) and not is_present(j, cur_ignored_cols):
                is_equivalent, is_empty = compute_intersection(n_row_long, cur_rows, &supports[j,0], row_buffer)

                if is_empty:
                    add_to_set(j, cur_ignored_cols)
                elif is_equivalent:
                    add_to_set(j, cur_cols)
                else:
                    NB_CHECKS += 1
                    idx_found = 0
                    found_one = False

                    for i in range(stack_length_new):
                        NB_SUB_CHECKS += 1
                        if is_subset(n_row_long, row_buffer, stack[i].rows):
                            found_one = True
                            break

                    if not found_one:
                        # canonic, add it to the stack
                        copy_set(n_row_long, row_buffer, stack[stack_length_new].rows)
                        stack[stack_length_new].current_column = j+1
                        stack_length_new += 1

        process()

        for i in range(stack_length, stack_length_new):
            copy_set(n_col_long, cur_cols, stack[i].cols)
            copy_set(n_col_long, cur_ignored_cols, stack[i].ignored_cols)
            add_to_set(stack[i].current_column - 1, stack[i].cols)

        stack_length = stack_length_new



    # todo = deque()
    # todo.append((set(range(0, matrix.shape[0])), set(), set(), {}, set(), 0))
    # while len(todo) != 0:
    #     rows, cols, P, N, parent, y = todo.pop()
    #     P = set(P)  # copy before modifying
    #     N = dict(N)
    #     todo_inside = []
    #     for j in range(y, matrix.shape[1]):
    #         if j not in cols and j not in P and j not in N:
    #             g = rows.intersection(supports[j])
    #             if len(g) == 0:
    #                 P.add(j)
    #             elif len(g) == len(rows):
    #                 cols.add(j)
    #             else:
    #                 NB_CHECKS += 1
    #                 idx_found = 0
    #                 found_one = False
    #
    #                 if not found_one:
    #                     for A, _, _, _, _, _ in todo:
    #                         NB_SUB_CHECKS += 1
    #                         if g.issubset(A):
    #                             found_one = True
    #                             break
    #                         idx_found += 1
    #
    #                 if not found_one:
    #                     for A, _ in todo_inside:
    #                         NB_SUB_CHECKS += 1
    #                         if g.issubset(A):
    #                             found_one = True
    #                             break
    #                         idx_found += 1
    #
    #                 if not found_one:
    #                     todo_inside.append((g, j))
    #                 else:
    #                     N[j] = idx_found
    #
    #     process(rows, cols)
    #     for g, j in todo_inside:
    #         h = cols.union({j})
    #         todo.append((g, h, P, {x:y for x, y in N.items() if y < len(todo)}, rows, j+1))

    for i in range(n_cols):
        free_stack_element(&stack[i])
    PyMem_Free(stack)
    PyMem_Free(row_buffer)
    PyMem_Free(cur_rows)
    PyMem_Free(cur_cols)
    PyMem_Free(cur_ignored_cols)

    return NB_CHECKS, NB_SUB_CHECKS