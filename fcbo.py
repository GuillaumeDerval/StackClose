from utils import compute_supports


def fcbo(matrix, process):

    NB_CANON = 0
    NB_INCL = 0
    NB_COLCHECK = 0

    n = matrix.shape[0]
    m = matrix.shape[1]

    supports = compute_supports(matrix)

    def canon_check(Nj, B, j):
        nonlocal NB_INCL
        NB_INCL += 1
        for i in range(0, j):
            if i in Nj and i not in B:
                return False
        return True

    def canon_check_eq(B, D, j):
        nonlocal NB_INCL
        NB_INCL += 1

        for i in range(0, j):
            if (i in B) != (i in D):
                return False
        return True

    def supported(A):
        nonlocal NB_CANON
        nonlocal NB_INCL
        NB_CANON += 1

        cols = set()
        for j in range(0, m):
            NB_INCL += 1
            if all(matrix[i, j] != 0 for i in A):
                cols.add(j)
        return cols

    def compute(A, B, y, N):
        nonlocal NB_COLCHECK

        if len(A) != 0:
            process(A, B)
        M = [s for s in N]
        todo = []
        for j in range(y, m):
            if j not in B:
                NB_COLCHECK += 1
                if canon_check(N[j], B, j):
                    C = A.intersection(supports[j])
                    D = supported(C)
                    if canon_check_eq(B, D, j):
                        todo.append((C, D, j))
                    else:
                        M[j] = D
        for C, D, j in todo:
            compute(C, D, j+1, M)

    compute(set(range(0, n)), supported(set(range(0, n))), 0, [set() for _ in range(0, m)])
    return NB_CANON, NB_INCL, NB_COLCHECK