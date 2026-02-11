import numpy as np

nodos = np.loadtxt('nodos.dat')
triangulos = np.loadtxt("triangulos.dat", dtype=int)

def punto_en_triangulo(p, p1, p2, p3):
    v1 = p2 - p1
    v2 = p3 - p1
    v = p - p1

    M = np.column_stack((v1, v2))
    try:
        alpha, beta = np.linalg.solve(M, v)
    except np.linalg.LinAlgError:
        return False

    return (alpha >= 0) and (beta >= 0) and (alpha + beta <= 1)

def encontrar_triangulo(p, nodos, triangulos):
    for i, tri in enumerate(triangulos):
        idx = tri - 1
        p1, p2, p3 = nodos[idx]

        if punto_en_triangulo(p, p1, p2, p3):
            return i, idx, (p1, p2, p3)

    return None, None, None

def funciones_forma(p1, p2, p3):
    A = np.array([
        [1, p1[0], p1[1]],
        [1, p2[0], p2[1]],
        [1, p3[0], p3[1]]
    ])

    N = []
    for i in range(3):
        b = np.zeros(3)
        b[i] = 1
        try:
            coef = np.linalg.solve(A, b)
        except np.linalg.LinAlgError:
            return None
        N.append(coef)

    return N

def f(x, y):
    return x**2 + y**2


def interpolar(p, nodos_triangulo):
    p1, p2, p3 = nodos_triangulo

    Z = np.array([
        f(p1[0], p1[1]),
        f(p2[0], p2[1]),
        f(p3[0], p3[1])
    ])

    N = funciones_forma(p1, p2, p3)
    if N is None:
        print("No se pudieron calcular las funciones de forma para el triángulo dado.")
        exit(1)

    x, y = p
    valor = 0.0
    for Zi, (a, b, c) in zip(Z, N):
        valor += Zi * (a + b*x + c*y)

    return valor

# Punto arbitrario
p = np.array([0.8, 0.2])

tri_id, idx, nodos_tri = encontrar_triangulo(p, nodos, triangulos)

if tri_id is None:
    print("El punto no está dentro del dominio")
else:
    z_interp = interpolar(p, nodos_tri)
    z_real = f(p[0], p[1])
    nodo_tri_encontrado = []
    for i, nodo in zip(idx, nodos_tri):
        nodo_tri_encontrado.append(int(i + 1))

    print(f"Triángulo encontrado: {tri_id + 1}")
    print(f"Nodos del triángulo: {nodo_tri_encontrado}")
    print("Interpolación:", z_interp)
    print("Valor exacto     :", z_real)
    print("Error            :", abs(z_interp - z_real))
