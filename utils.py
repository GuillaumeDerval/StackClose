# def compute_supports(matrix):
#     supports = [set() for _ in range(matrix.shape[1])]
#     for i in range(matrix.shape[0]):
#         for j in range(matrix.shape[1]):
#             if matrix[i,j]:
#                 supports[j].add(i)
#     return supports

def compute_supports(file):
    content = [[int(x) for x in line.strip().split(" ")] for line in open(file).read().strip().split("\n")]

    m = max(max(x) for x in content)+1
    supports = [set() for _ in range(m)]
    for idx, line in enumerate(content):
        for x in line:
            supports[x].add(idx)

    n_rows = len(content)
    n_cols = m
    return supports, n_rows, n_cols