import time
import numpy as np

from fcbo import fcbo
from inclose2_path import inclose2path
from inclose2_std import inclose2
from inclose4_path import inclose4path
from inclose4_std import inclose4
from inclose5_path import inclose5path
from inclose5_std import inclose5
from inclose5_path_cython import inclose5path as inclose5path_cython, inclose5 as inclose5_cython

def run_method_and_log(methodname, method, matrix):
    start = time.time()

    seen = []
    def process(*args):
        process.count += 1

    process.count = 0

    info_found = method(matrix, process)
    print("{:15} | {:10d} | {:8.2f} | {:12d} | {:12d} | {:12d}".format(
        methodname, process.count, time.time() - start, info_found[0], info_found[1], info_found[2]))

np.random.seed(16)
matrix = np.loadtxt("mushroom")
matrix = (np.random.randint(0, 5, (60, 60)) < 1).astype(np.bool).astype(np.double)

print(matrix.sum(), matrix.shape[0]*matrix.shape[1], matrix.sum()/(matrix.shape[0]*matrix.shape[1]))
print(matrix.shape)

print("{:15} | {:>10} | {:>8} | {:>12} | {:>12} | {:>12}".format("Method", "# Concepts", "Time (s)", "# canonicity", "# inclusion", "# colchecks"))

run_method_and_log("InClose5PathC", inclose5path_cython, matrix)
run_method_and_log("InClose5C", inclose5_cython, matrix)
print("----")
run_method_and_log("InClose5Path", inclose5path, matrix)
run_method_and_log("InClose5", inclose5, matrix)
print("----")
run_method_and_log("InClose4Path", inclose4path, matrix)
run_method_and_log("InClose4", inclose4, matrix)
print("----")
run_method_and_log("InClose2Path", inclose2path, matrix)
run_method_and_log("InClose2", inclose2, matrix)
print("----")
run_method_and_log("FCbO", fcbo, matrix)