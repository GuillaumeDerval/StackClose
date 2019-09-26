import sys
import numpy as np

matrix = np.loadtxt(sys.argv[1])

file = open(sys.argv[2], 'w')
for i in range(matrix.shape[0]):
    line = []
    for j in range(matrix.shape[1]):
        if matrix[i,j] != 0:
            line.append(j)

    file.write(" ".join(str(x) for x in line) + "\n")
file.close()