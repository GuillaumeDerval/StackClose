# distutils: language = c++
cimport cython
from set_utils cimport *
from cpython.mem cimport PyMem_Malloc, PyMem_Free

cdef struct stack_element:
    long* rows
    int row_set_size
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


def inclose5path(file, process, threshold=1):
    return __inclose5path(file, process, threshold)

def inclose5(file, process, threshold=1):
    return __inclose5(file, process, threshold)

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
cdef __inclose5path(file, process, int threshold):
    cdef unsigned long NB_CHECKS = 0
    cdef unsigned long NB_SUB_CHECKS = 0
    cdef unsigned long NB_COL_CHECKS = 0

    cdef int n_rows
    cdef int n_cols
    cdef long[:,:] supports
    cdef int n_row_long #= <int>ceil((<float>n_rows) / (8 * sizeof(long)))
    cdef int n_col_long #= <int>ceil((<float>n_cols) / (8 * sizeof(long)))

    supports, n_rows, n_cols, n_row_long, n_col_long = compute_supports(file)



    cdef stack_element* stack = <stack_element* >PyMem_Malloc(sizeof(stack_element) * n_rows)
    for i in range(n_cols):
        init_stack_element(&stack[i], n_row_long, n_col_long)

    cdef long* row_buffer = <long*>PyMem_Malloc(sizeof(long) * n_row_long)



    cdef int stack_length = 0

    # Add first element
    init_set(n_row_long, stack[0].rows)
    init_set(n_col_long, stack[0].cols)
    init_set(n_col_long, stack[0].ignored_cols)
    stack[0].current_column = 0
    stack[0].row_set_size = 0
    for i in range(n_rows):
        add_to_set(i, stack[0].rows)
        stack[0].row_set_size += 1
    stack_length += 1


    cdef long* cur_rows = <long*>PyMem_Malloc(sizeof(long) * n_row_long)
    cdef long* cur_cols = <long*>PyMem_Malloc(sizeof(long) * n_col_long)
    cdef long* cur_ignored_cols = <long*>PyMem_Malloc(sizeof(long) * n_col_long)

    cdef int cur_first_col
    cdef int cur_row_set_size
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
        cur_row_set_size = stack[stack_length].row_set_size
        stack_length_new = stack_length

        j = cur_first_col
        while j != n_cols:
            if not is_present(j, cur_ignored_cols):
                NB_COL_CHECKS += 1

                length = compute_intersection(n_row_long, cur_rows, &supports[j,0], row_buffer)

                if length == 0:
                    add_to_set(j, cur_ignored_cols)
                elif length == cur_row_set_size:
                    add_to_set(j, cur_cols)
                    add_to_set(j, cur_ignored_cols)
                else:
                    if length < threshold:
                        j += 1
                        continue

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
                        stack[stack_length_new].row_set_size = length
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
cdef __inclose5(file, process, int threshold):
    cdef unsigned long NB_CHECKS = 0
    cdef unsigned long NB_SUB_CHECKS = 0
    cdef unsigned long NB_COL_CHECKS = 0

    cdef int n_rows
    cdef int n_cols
    cdef long[:,:] supports
    cdef int n_row_long #= <int>ceil((<float>n_rows) / (8 * sizeof(long)))
    cdef int n_col_long #= <int>ceil((<float>n_cols) / (8 * sizeof(long)))

    supports, n_rows, n_cols, n_row_long, n_col_long = compute_supports(file)

    cdef stack_element* stack = <stack_element* >PyMem_Malloc(sizeof(stack_element) * n_cols)
    for i in range(n_cols):
        init_stack_element(&stack[i], n_row_long, n_col_long)

    cdef long* row_buffer = <long*>PyMem_Malloc(sizeof(long) * n_row_long)

    cdef int stack_length = 0

    # Add first element
    init_set(n_row_long, stack[0].rows)
    init_set(n_col_long, stack[0].cols)
    init_set(n_col_long, stack[0].ignored_cols)
    stack[0].current_column = 0
    stack[0].row_set_size = 0
    for i in range(n_rows):
        add_to_set(i, stack[0].rows)
        stack[0].row_set_size += 1
    stack_length += 1


    cdef long* cur_rows = <long*>PyMem_Malloc(sizeof(long) * n_row_long)
    cdef long* cur_cols = <long*>PyMem_Malloc(sizeof(long) * n_col_long)
    cdef long* cur_ignored_cols = <long*>PyMem_Malloc(sizeof(long) * n_col_long)

    cdef int cur_first_col
    cdef int cur_row_set_size
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
        cur_row_set_size = stack[stack_length].row_set_size
        stack_length_new = stack_length

        j = cur_first_col
        while j != n_cols:
            if not is_present(j, cur_ignored_cols):
                NB_COL_CHECKS += 1

                length = compute_intersection(n_row_long, cur_rows, &supports[j,0], row_buffer)

                if length == 0:
                    add_to_set(j, cur_ignored_cols)
                elif length == cur_row_set_size:

                    add_to_set(j, cur_cols)
                    add_to_set(j, cur_ignored_cols)
                else:
                    if length < threshold:
                        j += 1
                        continue

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
                        stack[stack_length_new].row_set_size = length
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
