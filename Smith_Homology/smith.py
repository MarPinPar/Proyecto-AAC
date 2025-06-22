from sympy import Matrix, eye, ZZ

# --- Funciones Auxiliares para Transformaciones de Fila y Columna ---

def multiplicar_izquierda_2(matriz, fila0_idx, fila1_idx, a, b, c, d):
    """
    Aplica una transformación a dos filas específicas de la matriz.
    Es como multiplicar la matriz por la izquierda por una matriz de 2x2.
    `matriz[fila0_idx]` se convierte en `a * matriz[fila0_idx] + b * matriz[fila1_idx]`
    `matriz[fila1_idx]` se convierte en `c * matriz[fila0_idx] + d * matriz[fila1_idx]`
    """
    # Itera sobre cada columna de las filas seleccionadas
    for j in range(matriz.cols):
        # Guarda los valores originales antes de modificarlos
        valor_fila0 = matriz[fila0_idx, j]
        valor_fila1 = matriz[fila1_idx, j]

        # Aplica la transformación lineal a las filas
        matriz[fila0_idx, j] = a * valor_fila0 + b * valor_fila1
        matriz[fila1_idx, j] = c * valor_fila0 + d * valor_fila1

def multiplicar_derecha_2(matriz, col0_idx, col1_idx, a, b, c, d):
    """
    Aplica una transformación a dos columnas específicas de la matriz.
    Es como multiplicar la matriz por la derecha por una matriz de 2x2.
    `matriz[:, col0_idx]` se convierte en `a * matriz[:, col0_idx] + c * matriz[:, col1_idx]`
    `matriz[:, col1_idx]` se convierte en `b * matriz[:, col0_idx] + d * matriz[:, col1_idx]`
    """
    # Itera sobre cada fila de las columnas seleccionadas
    for i in range(matriz.rows):
        # Guarda los valores originales antes de modificarlos
        valor_col0 = matriz[i, col0_idx]
        valor_col1 = matriz[i, col1_idx]

        # Aplica la transformación lineal a las columnas
        matriz[i, col0_idx] = a * valor_col0 + c * valor_col1
        matriz[i, col1_idx] = b * valor_col0 + d * valor_col1

# --- Función Principal: Forma Normal de Smith ---

def forma_normal_smith(matriz_entrada, dominio=ZZ):
    """
    Calcula la Forma Normal de Smith de una matriz.
    Retorna (S, D, T) donde S y T son matrices unimodulares y D es la forma de Smith,
    tal que S * Matriz_Original * T = D.
    El `dominio` por defecto es ZZ (números enteros).
    """
    # Convertir la entrada a un objeto Matrix de SymPy
    matriz = Matrix(matriz_entrada)

    # Inicializar las matrices de transformación S y T como matrices identidad
    # S transforma las filas (multiplicación por la izquierda)
    # T transforma las columnas (multiplicación por la derecha)
    s_transformacion = eye(matriz.rows)
    t_transformacion = eye(matriz.cols)

    # `ultima_col_procesada` ayuda a movernos por las columnas para no repetir trabajo
    ultima_col_procesada = -1

    # Bucle principal: Procesa la matriz para ponerla en forma de Smith
    # Itera a través de las filas de la matriz
    for i in range(matriz.rows):
        # Busca la primera columna no nula a partir de `ultima_col_procesada + 1`
        j_encontrada = -1 # Para rastrear la columna actual
        for j in range(ultima_col_procesada + 1, matriz.cols):
            if not matriz.col(j).is_zero: # Si la columna no es completamente cero
                j_encontrada = j
                break
        else: # Si no se encuentra ninguna columna no nula, la matriz está en forma de Smith (o es cero)
            break

        # Si el elemento pivote (matriz[i, j_encontrada]) es cero
        if matriz[i, j_encontrada] == 0:
            # Busca otra fila no nula en la misma columna
            fila_no_cero = -1
            for ii in range(matriz.rows):
                if matriz[ii, j_encontrada] != 0:
                    fila_no_cero = ii
                    break
            # Intercambia filas para llevar un elemento no nulo a la posición pivote
            multiplicar_izquierda_2(matriz, i, fila_no_cero, 0, 1, 1, 0)
            multiplicar_derecha_2(s_transformacion, i, fila_no_cero, 0, 1, 1, 0) # Actualiza S

        # Mueve la columna con el pivote a la posición actual `i`
        multiplicar_derecha_2(matriz, j_encontrada, i, 0, 1, 1, 0)
        multiplicar_izquierda_2(t_transformacion, j_encontrada, i, 0, 1, 1, 0) # Actualiza T
        # La columna activa ahora es `i`
        j_actual = i

        # Bucle interno: Reduce la fila y columna actuales a ceros (excepto el pivote)
        # y asegura la divisibilidad
        hubo_actualizacion = True
        while hubo_actualizacion:
            hubo_actualizacion = False

            # --- Reducción por Filas (Eliminación de elementos bajo el pivote) ---
            for fila_abajo in range(i + 1, matriz.rows):
                # Si el elemento es cero, no hay nada que hacer
                if matriz[fila_abajo, j_actual] == 0:
                    continue
                hubo_actualizacion = True # Se encontró un elemento no nulo, se necesita otra pasada

                # Si el elemento no es divisible por el pivote, usa el algoritmo extendido de Euclides
                # para hacer el pivote el MCD y luego cerar el otro elemento
                if dominio.rem(matriz[fila_abajo, j_actual], matriz[i, j_actual]) != 0:
                    coef1, coef2, g = dominio.gcdex(matriz[i, j_actual], matriz[fila_abajo, j_actual])
                    coef3 = dominio.quo(matriz[fila_abajo, j_actual], g)
                    coef4 = dominio.quo(matriz[i, j_actual], g)
                    multiplicar_izquierda_2(matriz, i, fila_abajo, coef1, coef2, -coef3, coef4)
                    multiplicar_derecha_2(s_transformacion, i, fila_abajo, coef4, -coef2, coef3, coef1) # Actualiza S

                # Finalmente, resta un múltiplo de la fila pivote para hacer el elemento cero
                coeficiente_division = dominio.quo(matriz[fila_abajo, j_actual], matriz[i, j_actual])
                multiplicar_izquierda_2(matriz, i, fila_abajo, 1, 0, -coeficiente_division, 1)
                multiplicar_derecha_2(s_transformacion, i, fila_abajo, 1, 0, coeficiente_division, 1) # Actualiza S

            # --- Reducción por Columnas (Eliminación de elementos a la derecha del pivote) ---
            for col_derecha in range(j_actual + 1, matriz.cols):
                # Si el elemento es cero, no hay nada que hacer
                if matriz[i, col_derecha] == 0:
                    continue
                hubo_actualizacion = True # Se encontró un elemento no nulo, se necesita otra pasada

                # Si el elemento no es divisible por el pivote, usa el algoritmo extendido de Euclides
                if dominio.rem(matriz[i, col_derecha], matriz[i, j_actual]) != 0:
                    coef1, coef2, g = dominio.gcdex(matriz[i, j_actual], matriz[i, col_derecha])
                    coef3 = dominio.quo(matriz[i, col_derecha], g)
                    coef4 = dominio.quo(matriz[i, j_actual], g)
                    multiplicar_derecha_2(matriz, j_actual, col_derecha, coef1, -coef3, coef2, coef4)
                    multiplicar_izquierda_2(t_transformacion, j_actual, col_derecha, coef4, coef3, -coef2, coef1) # Actualiza T

                # Finalmente, resta un múltiplo de la columna pivote para hacer el elemento cero
                coeficiente_division = dominio.quo(matriz[i, col_derecha], matriz[i, j_actual])
                multiplicar_derecha_2(matriz, j_actual, col_derecha, 1, -coeficiente_division, 0, 1)
                multiplicar_izquierda_2(t_transformacion, j_actual, col_derecha, 1, coeficiente_division, 0, 1) # Actualiza T

        # Actualiza la última columna procesada para la siguiente iteración del bucle principal
        ultima_col_procesada = j_actual

    # --- Post-procesamiento: Asegurar la propiedad de divisibilidad ---
    # Después de poner la matriz en forma diagonal, aseguramos que cada elemento divide al siguiente.
    # Esto se hace trabajando desde el final hacia el principio en la diagonal.
    for i1 in range(min(matriz.rows, matriz.cols)):
        for i0 in reversed(range(i1)):
            # Calcula el MCD y los coeficientes de Bezout
            coef1, coef2, g = dominio.gcdex(matriz[i0, i0], matriz[i1, i1])
            # Si el MCD es cero, no hay nada que hacer (uno o ambos elementos son cero)
            if g == 0:
                continue

            # Calcula coeficientes adicionales para las transformaciones
            coef3 = dominio.quo(matriz[i1, i1], g)
            coef4 = dominio.quo(matriz[i0, i0], g)

            # Aplica transformaciones de fila y columna para asegurar la divisibilidad
            # Estas transformaciones son más complejas, involucrando combinaciones
            # de filas/columnas para ajustar los elementos diagonales.
            multiplicar_izquierda_2(matriz, i0, i1, 1, coef2, coef3, coef2 * coef3 - 1)
            multiplicar_derecha_2(s_transformacion, i0, i1, 1 - coef2 * coef3, coef2, coef3, -1)
            multiplicar_derecha_2(matriz, i0, i1, coef1, 1 - coef1 * coef4, 1, -coef4)
            multiplicar_izquierda_2(t_transformacion, i0, i1, coef4, 1 - coef1 * coef4, 1, -coef1)

    # Retorna las matrices de transformación S, la matriz en forma de Smith (D) y la transformación T.
    return (s_transformacion, matriz, t_transformacion)

# --- Ejemplo de Uso ---

# Definición de la matriz de ejemplo
# Deberíamos añadir alguna otra y cambiar esa porque es la que viene en el código original.

mi_matriz = Matrix([
    [-1, -1, -1, 0, 0, 0, 0, 0, 0],
    [1, 0, 0, -1, -1, 0, 0, 0, 0],
    [0, 1, 0, 1, 0, -1, -1, 0, 0],
    [0, 0, 1, 0, 1, 1, 0, -1, 0],
    [0, 0, 0, 0, 0, 0, 1, 1, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, -1],
    [0, 0, 0, 0, 0, 0, 0, 0, 1]
])

# Calcular la Forma Normal de Smith
S, D, T = forma_normal_smith(mi_matriz)

# Imprimir los resultados
print("Mi Matriz Original:")
print(mi_matriz)
print("\nMatriz de Transformación S (para filas):")
print(S)
print("\nMatriz en Forma Normal de Smith (D):")
print(D)
print("\nMatriz de Transformación T (para columnas):")
print(T)

# Verificar la relación: S * Original * T = D
print("\nVerificación: S * Original * T = D")
print(S * mi_matriz * T)