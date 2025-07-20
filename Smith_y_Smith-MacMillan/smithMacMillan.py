# smith_macmillan.py
# Implementación del algoritmo iterativo para calcular la forma normal de Smith
# Basado en transformaciones de filas y columnas con operaciones enteras
# Esta versión busca coincidir con el comportamiento de forma_normal_smith de SymPy

import numpy as np
from math import gcd

def smith_macmillan(A_input):
    """
    Calcula la forma normal de Smith (Smith-MacMillan) para una matriz entera A.
    
    Parámetro:
        A_input (array-like): matriz de enteros (lista de listas, ndarray o similar)
    
    Retorna:
        ndarray: matriz diagonal D tal que D es equivalente a A bajo transformaciones
        elementales de filas y columnas, y satisface d_i | d_{i+1} en la diagonal.
    """
    A = np.array(A_input, dtype=int)
    A = A.copy()
    m, n = A.shape
    k = 0  # índice del pivote

    # Recorremos cada pivote desde (0,0) hasta min(m,n)
    while k < min(m, n):
        # Buscamos el menor valor absoluto distinto de cero en la submatriz A[k:, k:]
        sub = A[k:, k:]
        abs_sub = np.abs(sub)
        if not np.any(abs_sub):
            break  # Si el resto es nulo, terminamos

        # Encontramos el menor valor no nulo y lo llevamos a la posición (k,k)
        min_val = np.min(abs_sub[abs_sub > 0])
        i_min, j_min = np.argwhere(abs_sub == min_val)[0]
        i_min += k
        j_min += k

        A[[k, i_min]] = A[[i_min, k]]         # intercambiamos filas
        A[:, [k, j_min]] = A[:, [j_min, k]]   # intercambiamos columnas

        # Intentamos que A[k,k] divida a los demás elementos de su fila y columna
        done = False
        max_iter = 50
        iter_count = 0
        while not done and iter_count < max_iter:
            done = True
            for i in range(k + 1, m):
                if A[i, k] != 0:
                    g = gcd(abs(A[k, k]), abs(A[i, k]))
                    if g < abs(A[k, k]):
                        A[k] += A[i]
                        done = False
            for j in range(k + 1, n):
                if A[k, j] != 0:
                    g = gcd(abs(A[k, k]), abs(A[k, j]))
                    if g < abs(A[k, k]):
                        A[:, k] += A[:, j]
                        done = False
            iter_count += 1

        # Anulamos el resto de la columna y fila
        for i in range(m):
            if i != k and A[i, k] != 0:
                q = A[i, k] // A[k, k]
                A[i] -= q * A[k]
        for j in range(n):
            if j != k and A[k, j] != 0:
                q = A[k, j] // A[k, k]
                A[:, j] -= q * A[:, k]

        # Nos aseguramos de que el pivote sea positivo
        if A[k, k] < 0:
            A[k] = -A[k]

        k += 1  # pasamos al siguiente pivote

    # Extraemos la diagonal no nula y la ordenamos
    diag = [abs(A[i, i]) for i in range(min(m, n)) if A[i, i] != 0]
    diag.sort()

    # Forzar divisibilidad: d_i | d_{i+1}
    for i in range(len(diag) - 1):
        g = gcd(diag[i], diag[i + 1])
        diag[i + 1] = (diag[i + 1] * diag[i]) // g
        diag[i] = g

    # Post-procesamiento adicional para igualar a forma_normal_smith de SymPy
    # Si todos los elementos son 1 y no hay ninguna torsión significativa, descartamos
    if np.count_nonzero(A) == 0:
        diag = []


    # Construimos la matriz diagonal final
    D = np.zeros_like(A)
    for i in range(len(diag)):
        D[i, i] = diag[i]

    return D
