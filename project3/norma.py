import numpy as np
from numpy.linalg import inv

def condition_number(matrix_norm, M):
    if np.linalg.det(M) == 0:
        return float('inf')
    M_inv = inv(M)
    return matrix_norm(M) * matrix_norm(M_inv)

def matrix_norm_1(M):
    return np.max(np.sum(np.abs(M), axis=0))



def matrix_norm_2(M):
    matrix = M.T @ M
    eigenvalues = np.linalg.eigvals(matrix)
    return np.sqrt(max(eigenvalues))


def matrix_norm_p(M, p):
    return (sum(abs(M[i][j])**p for i in range(len(M)) for j in range(len(M[0])))) ** (1/p)


def matrix_norm_inf(M):
    return max(np.sum(np.abs((M)), axis=1))



M = np.random.randint(-10, 10, size=(2, 2))
print("\nMacierz M=\n", M)
print("zaimplementowana norma 1:", matrix_norm_1(M))
print("numpy norma 1:", np.linalg.norm(M, 1))
print("zaimplementowana norma 2:", matrix_norm_2(M))
print("numpy norma 2:", np.linalg.norm(M, 2))
print("zaimplementowana norma inf:", matrix_norm_inf(M))
print("numpy norma inf:", np.linalg.norm(M, np.inf))
p = 2
print("zaimplementowana norma p:", matrix_norm_p(M, p))
print("numpy norma p :", np.linalg.norm(M, 'fro'))
print("współczynnik uwarunkowania:", condition_number(matrix_norm_2, M))
print("numpy współczynnik uwarunkowania:", np.linalg.cond(M, 2))