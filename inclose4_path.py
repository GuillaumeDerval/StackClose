from collections import deque

from utils import compute_supports


def inclose4path(matrix, process):
    NB_CHECKS = 0
    NB_SUB_CHECKS = 0
    NB_COL_CHECKS = 0

    supports = compute_supports(matrix)

    max_todo_size = 0
    todo = deque()
    todo.append((set(range(0, matrix.shape[0])), set(), set(), 0))
    while len(todo) != 0:
        max_todo_size = max(max_todo_size, len(todo))
        rows, cols, P, y = todo.pop()
        P = set(P)  # copy before modifying
        todo_inside = []
        for j in range(y, matrix.shape[1]):
            if j not in cols and j not in P:
                NB_COL_CHECKS += 1
                g = rows.intersection(supports[j])
                if len(g) == 0:
                    P.add(j)
                elif len(g) == len(rows):
                    cols.add(j)
                else:
                    NB_CHECKS += 1
                    found_one = False

                    if not found_one:
                        for A, _, _, _ in todo:
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

        process(rows, cols)
        for g, j in todo_inside:
            h = cols.union({j})
            todo.append((g, h, P, j+1))

    return NB_CHECKS, NB_SUB_CHECKS, NB_COL_CHECKS