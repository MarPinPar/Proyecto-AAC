import numpy as np
from sympy import Matrix
from smith import forma_normal_smith
from smithMacMillan import smith_macmillan

def extraer_diagonal_efectiva(mat):
    d = [abs(mat[i, i]) for i in range(min(mat.shape)) if mat[i, i] != 0]
    d.sort()
    return d

def es_forma_smith_valida(diag):
    for i in range(len(diag) - 1):
        if diag[i] == 0:
            continue
        if diag[i+1] % diag[i] != 0:
            return False
    return True

# Convierte elementos numpy a int normales
def limpiar_lista_numeros(lista):
    return [int(x) for x in lista]

def comparar(matriz_sym, nombre):
    print(f"\n=== {nombre} ===")
    matriz_np = np.array(matriz_sym.tolist(), dtype=int)

    resultado_np = smith_macmillan(matriz_np)
    diag_np = extraer_diagonal_efectiva(resultado_np)
    print("smith_macmillan:", limpiar_lista_numeros(diag_np))

    _, resultado_sym, _ = forma_normal_smith(matriz_sym)
    diag_sym = extraer_diagonal_efectiva(resultado_sym)
    print("forma_normal_smith:", limpiar_lista_numeros(diag_sym))

    if diag_np == diag_sym:
        print("Coinciden los invariantes elementales.")
    else:
        print("Diferencias encontradas en los invariantes.")

    if not es_forma_smith_valida(diag_np):
        print("La diagonal de smith_macmillan no cumple d_i | d_{i+1}")
    if not es_forma_smith_valida(diag_sym):
        print("La diagonal de forma_normal_smith no cumple d_i | d_{i+1}")

if __name__ == "__main__":
    # Definiciones locales de las matrices de prueba
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
        comparar(m, f"m_{i}")
