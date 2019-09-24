from collections import deque
from typing import Set

from utils import compute_supports

def inclose2(matrix, process):

    def is_canonical(supports, g: Set[int], cols: Set[int], j):
        is_canonical.n_checks += 1
        for jj in range(0, j):
            if jj not in cols:
                is_canonical.n_sub_checks += 1
                if g.issubset(supports[jj]):
                    return False
        return True

    is_canonical.n_checks = 0
    is_canonical.n_sub_checks = 0
    NB_COL_CHECKS = 0

    supports = compute_supports(matrix)

    max_todo_size = 0
    todo = deque()
    todo.append((set(range(0, matrix.shape[0])), set(), 0))
    while len(todo) != 0:
        max_todo_size = max(max_todo_size, len(todo))
        rows, cols, y = todo.pop()
        todo_inside = []
        for j in range(y, matrix.shape[1]):
            if j not in cols:
                NB_COL_CHECKS += 1
                g = rows.intersection(supports[j])
                if len(g) == len(rows):
                    cols.add(j)
                else:
                    canonical = is_canonical(supports, g, cols, j)
                    if canonical:
                        todo_inside.append((g, j))

        if len(rows) != 0:
            process(rows, cols)
        #todo_inside = reversed(todo_inside)
        for g, j in todo_inside:
            h = cols.union({j})
            todo.append((g, h, j+1))

    return is_canonical.n_checks, is_canonical.n_sub_checks, NB_COL_CHECKS