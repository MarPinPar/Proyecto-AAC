from sympy import Matrix
from smith import forma_normal_smith

def probar_caso(matriz):
    print("\nProbando matriz:")
    print(matriz)
    S, D, T = forma_normal_smith(matriz)
    print("\nForma Normal de Smith D:")
    print(D)
    verificacion = S * matriz * T
    assert verificacion.equals(D), "Error: La verificación S * A * T = D falló"
    print("Verificación exitosa: S * A * T = D")

if __name__ == "__main__":
    # Caso base que ya tienes
    matriz1 = Matrix([
        [-1, -1, -1, 0, 0, 0, 0, 0, 0],
        [1,  0,  0, -1, -1, 0, 0, 0, 0],
        [0,  1,  0,  1,  0, -1, -1, 0, 0],
        [0,  0,  1,  0,  1,  1,  0, -1, 0],
        [0,  0,  0,  0,  0,  0,  1, 1, 0],
        [0,  0,  0,  0,  0,  0,  0, 0, -1],
        [0,  0,  0,  0,  0,  0,  0, 0, 1]
    ])
    probar_caso(matriz1)

    # Aquí agrega tus matrices para probar:
    matriz2 = Matrix([
        [2, 4],
        [6, 8]
    ])
    probar_caso(matriz2)

    matriz3 = Matrix([
        [0, 1, 2],
        [3, 4, 5],
        [6, 7, 8]
    ])
    probar_caso(matriz3)

    # Añade más según tus casos
