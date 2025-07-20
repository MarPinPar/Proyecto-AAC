import numpy as np
from sympy import Matrix
from smithMacMillan import smith_macmillan  # asegúrate que el archivo se llame así

def es_forma_normal_de_smith(diag):
    # Verifica que la lista diag cumpla d_i | d_{i+1}
    for i in range(len(diag) - 1):
        if diag[i] == 0:
            continue
        if diag[i+1] % diag[i] != 0:
            return False
    return True

def probar(matriz_sym, nombre):
    print(f"\n=== {nombre} ===")
    matriz_np = np.array(matriz_sym.tolist(), dtype=int)
    print("Matriz original:")
    print(matriz_np)

    resultado = smith_macmillan(matriz_np)
    print("Resultado Smith-MacMillan:")
    print(resultado)

    # Extraer y verificar la diagonal
    min_dim = min(resultado.shape)
    diag = [abs(resultado[i, i]) for i in range(min_dim) if resultado[i, i] != 0]
    diag_limpia = [int(x) for x in diag]
    print("Diagonal no nula:", diag_limpia)

    if es_forma_normal_de_smith(diag):
        print("Es una forma normal de Smith válida.")
    else:
        print("No cumple d_i | d_{i+1}.")

if __name__ == "__main__":
    # Matrices de prueba
    m_1 = Matrix([
        [6, 2, 3],
        [2, 4, 0],
        [3, 0, 1]
    ])

    m_2 = Matrix([
        [0, -1, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, -1],
        [-1, 0, 1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0],
        [0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0],
        [1, 1, 0, 1, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, -1, 1, 1, 0, 0, 1, 0],
        [0, 0, 0, 0, 0, 1, 1, 0, -1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 1, 1, 0, 1],
        [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0]
    ])

    m_3 = Matrix([
        [1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1],
        [-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0],
        [0, -1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
        [0, -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0],
        [0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1],
        [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 1],
        [0, 0, 0, -1, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0],
        [0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
        [0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0, 0, 0, -1, 0],
        [0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0, -1, 0],
        [0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0]
    ])

    m_4 = Matrix([
        [0, -1, -1, 0, -1, 0, -1, 0, -1, 0, 0, 0],
        [-1, 0, 1, 0, 0, 0, 0, 0, 0, -1, 0, -1],
        [1, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, 1, -1, 0, 0, 0, 0, -1, 1],
        [0, 0, 0, 0, 0, 1, 1, -1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0]
    ])

    m_5 = Matrix([
        [1, 0, 0, 0, 0, 0],
        [-1, 1, 0, 0, 0, 0],
        [1, 0, 0, 0, -1, 0],
        [0, 1, 0, 0, 0, 0],
        [0, -1, 1, 0, 0, 0],
        [0, 0, 1, 0, 0, 0],
        [0, 0, -1, 1, 0, 0],
        [0, 0, 0, 1, 0, 0],
        [0, 0, 0, 1, 0, 0],
        [0, 0, 0, -1, 1, 0],
        [0, 0, 0, 0, -1, -1],
        [0, 0, 0, 0, 0, 1],
        [0, 0, 0, 0, 0, 1]
    ])

    matrices = [m_1, m_2, m_3, m_4, m_5]
    for i, m in enumerate(matrices, 1):
        probar(m, f"m_{i}")
