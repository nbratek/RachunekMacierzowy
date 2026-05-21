import numpy as np
import matplotlib.pyplot as plt

def power_method(A, p=2, epsilon=1e-4, max_iter=10000):
    n = A.shape[0]
    while True:
        z = np.random.uniform(1e-6, 1-1e-6, n)
        w = A @ z
        eigenvalue = w[np.argmax(np.abs(w))]
        if np.linalg.norm(w - eigenvalue * z, ord=p) >= 1e-8:
            break
    errors = []
    for i in range(max_iter):
        w = A @ z
        eigenvalue = w[np.argmax(np.abs(w))]
        error = np.linalg.norm(w - eigenvalue * z, ord=p)
        errors.append(error)
        if error < epsilon:
            break
        z = w / eigenvalue
    v = z / np.linalg.norm(z, 2)
    return eigenvalue, v, errors


A = np.random.rand(3, 3)
print("A =\n", A)

AAT = A @ A.T
eigvals, U = np.linalg.eigh(AAT)
D = np.diag(np.sqrt(eigvals[::-1]))
U = U[:, ::-1]
D_inv = np.linalg.inv(D)
V = A.T @ U @ D_inv

print("U =\n", U)
print("D =\n", D)
print("V =\n", V)

A_def = A.copy()
for k in range(3):
    eigenvalue_def, v_def, _ = power_method(A_def, p=2, epsilon=1e-4)
    for p in [1, 2, 3, 4, np.inf]:
        eigenvalue, v, errors = power_method(A_def, p=p, epsilon=1e-4)
        plt.plot(errors, marker='o', markersize=3)
        plt.yscale('log')
        plt.xlabel('Iteracja')
        plt.ylabel('Błąd')
        plt.title(f'Zbieżność metody potęgowej, λ_{k + 1} ≈ {eigenvalue:.3f}, norma p={p}')
        plt.grid(True)
        plt.savefig(f'zbieznosc_{k + 1}_p{p}.png')
        plt.clf()
    A_def = A_def - eigenvalue_def * np.outer(v_def, v_def)




U_numpy, S_numpy, VT_numpy = np.linalg.svd(A)
A_numpy = U_numpy @ np.diag(S_numpy) @ VT_numpy
diff = (U @ D @ V.T - A_numpy).flatten()
for p in [1, 2, 3, 4, np.inf]:
    print(f"p={p}: {np.linalg.norm(diff, ord=p):.3e}")