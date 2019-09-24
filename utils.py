def compute_supports(matrix):
    supports = [set() for _ in range(matrix.shape[1])]
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if matrix[i,j]:
                supports[j].add(i)
    return supports