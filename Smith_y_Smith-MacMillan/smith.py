from sympy import Matrix, eye, ZZ

def swap_filas(matriz, i, j):
    matriz[i, :], matriz[j, :] = matriz[j, :], matriz[i, :]

def swap_columnas(matriz, i, j):
    matriz[:, i], matriz[:, j] = matriz[:, j], matriz[:, i]

def multiplicar_izquierda_2(matriz, fila0_idx, fila1_idx, a, b, c, d):
    temp0 = []
    temp1 = []
    for j in range(matriz.cols):
        val0 = matriz[fila0_idx, j]
        val1 = matriz[fila1_idx, j]
        temp0.append(a * val0 + b * val1)
        temp1.append(c * val0 + d * val1)
    for j in range(matriz.cols):
        matriz[fila0_idx, j] = temp0[j]
        matriz[fila1_idx, j] = temp1[j]

def multiplicar_derecha_2(matriz, col0_idx, col1_idx, a, b, c, d):
    temp0 = []
    temp1 = []
    for i in range(matriz.rows):
        val0 = matriz[i, col0_idx]
        val1 = matriz[i, col1_idx]
        temp0.append(a * val0 + b * val1)
        temp1.append(c * val0 + d * val1)
    for i in range(matriz.rows):
        matriz[i, col0_idx] = temp0[i]
        matriz[i, col1_idx] = temp1[i]

def division_segura(a, b, dominio):
    return dominio.quo(a, b) if b != 0 else 0

def forma_normal_smith(matriz_entrada, dominio=ZZ, verbose=False):
    if verbose:
        print("Iniciando cálculo de Forma Normal de Smith...")
    matriz = Matrix(matriz_entrada).copy()
    s_transformacion = eye(matriz.rows)
    t_transformacion = eye(matriz.cols)
    ultima_col_procesada = -1

    for i in range(matriz.rows):
        if verbose:
            print(f"\nProcesando fila {i}...")
        j_encontrada = -1
        for j in range(ultima_col_procesada + 1, matriz.cols):
            if not matriz.col(j).is_zero:
                j_encontrada = j
                break
        else:
            if verbose:
                print("  No se encontró columna no nula. Terminando.")
            break

        if verbose:
            print(f"  Columna pivote encontrada: {j_encontrada}")

        if matriz[i, j_encontrada] == 0:
            fila_no_cero = -1
            for ii in range(matriz.rows):
                if matriz[ii, j_encontrada] != 0:
                    fila_no_cero = ii
                    break
            if verbose:
                print(f"  Intercambiando fila {i} con fila {fila_no_cero} para obtener pivote no nulo")
            swap_filas(matriz, i, fila_no_cero)
            swap_filas(s_transformacion, i, fila_no_cero)

        swap_columnas(matriz, j_encontrada, i)
        swap_columnas(t_transformacion, j_encontrada, i)
        j_actual = i
        if verbose:
            print(f"  Pivote actual en ({i}, {j_actual}) = {matriz[i, j_actual]}")

        filas_reducidas = set()
        columnas_reducidas = set()

        hubo_actualizacion = True
        while hubo_actualizacion:
            if verbose:
                print("  Iteración de reducción en curso...")
            hubo_actualizacion = False

            for fila_abajo in range(i + 1, matriz.rows):
                if fila_abajo in filas_reducidas:
                    continue
                valor_antes = matriz[fila_abajo, j_actual]
                if valor_antes == 0:
                    filas_reducidas.add(fila_abajo)
                    continue
                if verbose:
                    print(f"    Reduciendo fila {fila_abajo} usando fila {i}")
                hubo_actualizacion = True

                if matriz[i, j_actual] != 0 and dominio.rem(valor_antes, matriz[i, j_actual]) != 0:
                    if verbose:
                        print("      División no exacta: aplicando gcdex por filas")
                    coef1, coef2, g = dominio.gcdex(matriz[i, j_actual], valor_antes)
                    coef3 = division_segura(valor_antes, g, dominio)
                    coef4 = division_segura(matriz[i, j_actual], g, dominio)
                    multiplicar_izquierda_2(matriz, i, fila_abajo, coef1, coef2, -coef3, coef4)
                    multiplicar_izquierda_2(s_transformacion, i, fila_abajo, coef1, coef2, -coef3, coef4)

                if matriz[i, j_actual] != 0:
                    coef_div = division_segura(matriz[fila_abajo, j_actual], matriz[i, j_actual], dominio)
                    multiplicar_izquierda_2(matriz, i, fila_abajo, 1, 0, -coef_div, 1)
                    multiplicar_izquierda_2(s_transformacion, i, fila_abajo, 1, 0, -coef_div, 1)

                if matriz[fila_abajo, j_actual] == 0:
                    filas_reducidas.add(fila_abajo)

            for col_derecha in range(j_actual + 1, matriz.cols):
                if col_derecha in columnas_reducidas:
                    continue
                valor_antes = matriz[i, col_derecha]
                if valor_antes == 0 or matriz[i, j_actual] == 0:
                    columnas_reducidas.add(col_derecha)
                    continue
                if verbose:
                    print(f"    Reduciendo columna {col_derecha} usando columna {j_actual}")
                hubo_actualizacion = True

                if dominio.rem(valor_antes, matriz[i, j_actual]) != 0:
                    if verbose:
                        print("      División no exacta: aplicando gcdex por columnas")
                    coef1, coef2, g = dominio.gcdex(matriz[i, j_actual], valor_antes)
                    coef3 = division_segura(valor_antes, g, dominio)
                    coef4 = division_segura(matriz[i, j_actual], g, dominio)
                    multiplicar_derecha_2(matriz, j_actual, col_derecha, coef1, coef2, -coef3, coef4)
                    multiplicar_derecha_2(t_transformacion, j_actual, col_derecha, coef1, coef2, -coef3, coef4)

                coef_div = division_segura(matriz[i, col_derecha], matriz[i, j_actual], dominio)
                multiplicar_derecha_2(matriz, j_actual, col_derecha, 1, -coef_div, 0, 1)
                multiplicar_derecha_2(t_transformacion, j_actual, col_derecha, 1, -coef_div, 0, 1)

                if matriz[i, col_derecha] == 0:
                    columnas_reducidas.add(col_derecha)

        ultima_col_procesada = j_actual

    if verbose:
        print("Cálculo completado.")
    return (s_transformacion, matriz, t_transformacion)

def normalizar_forma_smith(D):
    from math import gcd
    from functools import reduce

    # Extraer diagonal significativa
    diag = [abs(D[i, i]) for i in range(min(D.rows, D.cols)) if D[i, i] != 0]
    diag.sort()

    # Ajustar para que d_i | d_{i+1}
    for i in range(len(diag) - 1):
        g = gcd(diag[i], diag[i+1])
        if g != diag[i]:
            diag[i+1] = (diag[i+1] * diag[i]) // g
            diag[i] = g

    # Reconstruir matriz diagonal
    D_final = Matrix.zeros(D.rows, D.cols)
    for i, val in enumerate(diag):
        D_final[i, i] = val

    return D_final
