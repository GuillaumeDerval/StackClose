from collections import deque

from utils import compute_supports


def inclose2path(file, process, thresold=1):
    NB_CHECKS = 0
    NB_SUB_CHECKS = 0
    NB_COL_CHECKS = 0

    supports, n_rows, n_cols = compute_supports(file)

    max_todo_size = 0
    todo = deque()
    todo.append((set(range(0, n_rows)), set(), 0))
    while len(todo) != 0:
        max_todo_size = max(max_todo_size, len(todo))
        rows, cols, y = todo.pop()
        todo_inside = []
        for j in range(y, n_cols):
            if j not in cols:
                NB_COL_CHECKS += 1
                g = rows.intersection(supports[j])
                if len(g) == len(rows):
                    cols.add(j)
                else:
                    if len(g) < thresold:
                        continue

                    NB_CHECKS += 1
                    found_one = False

                    if not found_one:
                        for A, _, jj in todo:
                            NB_SUB_CHECKS += 1
                            if g.issubset(A):
                                found_one = True
                                break

                    if not found_one:
                        for A, _ in todo_inside:
                            NB_SUB_CHECKS += 1
                            if g.issubset(A):
                                found_one = True
                                break

                    if not found_one:
                        todo_inside.append((g, j))

        if len(rows) != 0:
            process(rows, cols)
        for g, j in todo_inside:
            h = cols.union({j})
            todo.append((g, h, j+1))

    return NB_CHECKS, NB_SUB_CHECKS, NB_COL_CHECKS