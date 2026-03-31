import numpy as np

def lu_no_pivot(A):
    n = A.shape[0]
    A = A.astype(float)

    L = np.eye(n)
    U = A.copy()

    for k in range(n):
        for i in range(k + 1, n):
            factor = U[i, k] / U[k, k]
            L[i, k] = factor
            U[i, :] -= factor * U[k, :]

    return L, U

def lu_partial_pivot(A):
    n = A.shape[0]
    A = A.astype(float)

    P = np.eye(n)
    L = np.eye(n)
    U = A.copy()

    for k in range(n):
        max_row = np.argmax(abs(U[k:, k])) + k

        if max_row != k:
            U[[k, max_row]] = U[[max_row, k]]
            P[[k, max_row]] = P[[max_row, k]]

            if k > 0:
                L[[k, max_row], :k] = L[[max_row, k], :k]

        for i in range(k + 1, n):
            factor = U[i, k] / U[k, k]
            L[i, k] = factor
            U[i, :] -= factor * U[k, :]

    return P, L, U