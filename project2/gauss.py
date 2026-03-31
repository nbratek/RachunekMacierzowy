import numpy as np

def gauss_no_pivot_unit_diag(A, b):
    n = A.shape[0]
    A = A.astype(float)
    b = b.astype(float)

    for k in range(n):
        pivot = A[k, k]
        
        A[k, :] = A[k, :] / pivot
        b[k] = b[k] / pivot

        for i in range(k + 1, n):
            factor = A[i, k]
            A[i, :] -= factor * A[k, :]
            b[i] -= factor * b[k]

    x = np.zeros(n)
    for i in reversed(range(n)):
        x[i] = b[i] - np.dot(A[i, i+1:], x[i+1:])

    return x


def gauss_partial_pivot(A, b):
    n = A.shape[0]
    A = A.astype(float)
    b = b.astype(float)

    for k in range(n):
        max_row = np.argmax(abs(A[k:, k])) + k

        if max_row != k:
            A[[k, max_row]] = A[[max_row, k]]
            b[[k, max_row]] = b[[max_row, k]]

        for i in range(k + 1, n):
            factor = A[i, k] / A[k, k]
            A[i, :] -= factor * A[k, :]
            b[i] -= factor * b[k]

    x = np.zeros(n)
    for i in reversed(range(n)):
        x[i] = (b[i] - np.dot(A[i, i+1:], x[i+1:])) / A[i, i]

    return x