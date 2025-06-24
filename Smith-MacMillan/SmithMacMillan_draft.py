'''
Implementación del algoritmo de Smith-MacMillan para matrices enteras.
Código en proceso, es posible que no funcione correctamente en algunos casos de prueba.
'''


import numpy as np
import math

def smith_macmillan(A_input):
    # A es una matriz de tamaño m x n con entradas enteras
    A = np.array(A_input, dtype=int) # Convertir a numpy array para operaciones eficientes

    # Si A es la matriz nula:
    # Si A es la matriz nula, retornar A.
    if np.all(A == 0):
        return A

    m, n = A.shape

    # Paso 1: Llevar el menor valor absoluto no nulo a la posición (1,1)
    # Encontrar el elemento no nulo de menor valor absoluto
    min_val = float('inf')
    min_r, min_c = -1, -1
    for r in range(m):
        for c in range(n):
            if A[r, c] != 0 and abs(A[r, c]) < min_val:
                min_val = abs(A[r, c])
                min_r, min_c = r, c

    # Usar permutaciones de filas y columnas para moverlo a A[0][0] (índice 0 en Python)
    if min_r != 0:
        A[[0, min_r]] = A[[min_r, 0]] # Intercambiar fila 0 y min_r
    if min_c != 0:
        A[:, [0, min_c]] = A[:, [min_c, 0]] # Intercambiar columna 0 y min_c

    # Mientras A[0][0] no divida a todos los elementos de la primera fila y columna:
    # Repetir mientras A[0,0] no divida completamente la primera fila o columna.
    while True:
        changed = False

        # Para cada elemento A[0][j] (j > 0):
        # Aplicar algoritmo de Euclides para reducir elementos de la primera fila.
        for j in range(1, n):
            if A[0, 0] != 0 and A[0, j] % A[0, 0] != 0:
                q = A[0, j] // A[0, 0]
                A[0, j] -= q * A[0, 0]
                changed = True

        # Para cada elemento A[i][0] (i > 0):
        # Aplicar algoritmo de Euclides para reducir elementos de la primera columna.
        for i in range(1, m):
            if A[0, 0] != 0 and A[i, 0] % A[0, 0] != 0:
                q = A[i, 0] // A[0, 0]
                A[i, 0] -= q * A[0, 0]
                changed = True
        
        # Después de las operaciones, si A[0,0] se ha convertido en 0, 
        # necesitamos reubicar el menor elemento no nulo en A[0,0]
        if A[0,0] == 0 and not np.all(A == 0):
            # Encontrar el elemento no nulo de menor valor absoluto
            min_val = float('inf')
            min_r, min_c = -1, -1
            for r_inner in range(m):
                for c_inner in range(n):
                    if A[r_inner, c_inner] != 0 and abs(A[r_inner, c_inner]) < min_val:
                        min_val = abs(A[r_inner, c_inner])
                        min_r, min_c = r_inner, c_inner

            # Si se encontró un elemento no nulo, moverlo a A[0][0] y continuar el bucle
            if min_r != -1:
                if min_r != 0:
                    A[[0, min_r]] = A[[min_r, 0]]
                if min_c != 0:
                    A[:, [0, min_c]] = A[:, [min_c, 0]]
                changed = True
            else: # Si todos los elementos son 0, la matriz es nula
                return np.zeros_like(A)

        if not changed:
            break

    # Ahora A[0][0] divide todos los elementos de la primera fila y columna
    # Paso 3: Asegurar que A[0][0] divide al resto de elementos de la submatriz
    # Iterar sobre la submatriz y ajustar si un elemento no es divisible por A[0,0].
    if A[0,0] != 0:
        for i in range(1, m):
            for j in range(1, n):
                # Mientras A[i][j] no sea divisible por A[0][0]:
                # Sumar fila 0 a fila i o columna 0 a columna j.
                while A[i, j] % A[0, 0] != 0:
                    # Estrategia: sumar la fila 0 a la fila i
                    # o la columna 0 a la columna j, lo que sea más sencillo.
                    # Aquí sumamos la fila 0 a la fila i para cambiar A[i,j]
                    A[i, :] += A[0, :]
    
    # Paso 4: A[0][0] divide todos los elementos, aplicar operaciones de fila/columna
    # Hacer que los elementos de la primera fila y columna (excepto A[0,0]) sean cero.
    for i in range(1, m):
        # A[i][0] - (A[i][0] / A[0][0]) * A[0][0] = 0 si es divisible
        if A[0,0] != 0:
            A[i, :] -= (A[i, 0] // A[0, 0]) * A[0, :]
    for j in range(1, n):
        # A[0][j] - (A[0][j] / A[0][0]) * A[0][0] = 0 si es divisible
        if A[0,0] != 0:
            A[:, j] -= (A[0, j] // A[0, 0]) * A[:, 0]

    # Paso 5: Separar la matriz en la forma A = diag(d1, 0) o A'
    # Extraer el elemento pivote d1.
    d1 = A[0, 0]
    # Extraer la submatriz A' (eliminando la primera fila y columna).
    A_prime = A[1:, 1:]

    # Paso 6: Aplicar el algoritmo recursivamente a A'
    # Llamada recursiva a SmithMacMillan para la submatriz.
    A_prime_smith = smith_macmillan(A_prime)

    # Construir matriz diagonal completa
    # Construir la matriz de resultado final.
    result_matrix = np.zeros_like(A)
    result_matrix[0, 0] = d1
    result_matrix[1:, 1:] = A_prime_smith

    return result_matrix