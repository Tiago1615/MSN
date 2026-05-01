import numpy as np


# ============================================================
# CONFIGURACIÓN DEL PROBLEMA
# ============================================================

L = 10.0
Nelements = 100
c = 1.0

dirichlet_conditions = {
    "left": 0.0
}

neumann_conditions = {
    "right": 1.0
}


def f(x):
    return 2.0 * x


# ============================================================
# MALLA
# ============================================================

def crear_malla(L, Nelements):
    Nnodes = Nelements + 1
    x = np.linspace(0.0, L, Nnodes)
    return x, Nnodes


# ============================================================
# CUADRATURA Y FUNCIONES DE FORMA
# ============================================================

def gauss_2p():
    xi = np.array([
        -1.0 / np.sqrt(3.0),
         1.0 / np.sqrt(3.0)
    ])

    w = np.array([1.0, 1.0])

    return xi, w


def funciones_forma_lineales(xi):
    N = np.array([
        (1.0 - xi) / 2.0,
        (1.0 + xi) / 2.0
    ])

    dN_dxi = np.array([
        -0.5,
         0.5
    ])

    return N, dN_dxi


# ============================================================
# MATRICES LOCALES
# ============================================================

def calcular_matriz_local(x1, x2, c, f):
    h = x2 - x1
    J = h / 2.0

    Ae = np.zeros((2, 2))
    Be = np.zeros(2)

    xi_gauss, pesos = gauss_2p()

    for xi, w in zip(xi_gauss, pesos):

        N, dN_dxi = funciones_forma_lineales(xi)

        x_phys = (x1 + x2) / 2.0 + J * xi

        dN_dx = dN_dxi / J

        for a in range(2):

            for b in range(2):

                # Término de rigidez
                Ae[a, b] += w * dN_dx[a] * dN_dx[b] * J

                # Término de reacción
                Ae[a, b] += w * c * N[a] * N[b] * J

            # Vector de cargas
            Be[a] += w * f(x_phys) * N[a] * J

    return Ae, Be


# ============================================================
# ENSAMBLAJE GLOBAL
# ============================================================

def ensamblar_sistema(x, Nelements, c, f):
    Nnodes = len(x)

    A = np.zeros((Nnodes, Nnodes))
    B = np.zeros(Nnodes)

    for e in range(Nelements):

        x1 = x[e]
        x2 = x[e + 1]

        Ae, Be = calcular_matriz_local(x1, x2, c, f)

        nodos = [e, e + 1]

        for a in range(2):

            I = nodos[a]
            B[I] += Be[a]

            for b in range(2):

                Jg = nodos[b]
                A[I, Jg] += Ae[a, b]

    return A, B


# ============================================================
# CONDICIONES DE CONTORNO
# ============================================================

def obtener_nodo_lado(lado, Nnodes):
    if lado == "left":
        return 0

    if lado == "right":
        return Nnodes - 1

    raise ValueError(f"Lado no reconocido: {lado}")


def aplicar_neumann(B, neumann_conditions, Nnodes):
    for lado, valor in neumann_conditions.items():

        nodo = obtener_nodo_lado(lado, Nnodes)

        if lado == "right":
            B[nodo] += valor

        elif lado == "left":
            B[nodo] -= valor

    return B


def aplicar_dirichlet(A, B, dirichlet_conditions, Nnodes):
    nodos_dirichlet = {}

    for lado, valor in dirichlet_conditions.items():
        nodo = obtener_nodo_lado(lado, Nnodes)
        nodos_dirichlet[nodo] = valor

    # Corrección del segundo miembro
    for nodo, valor in nodos_dirichlet.items():
        B -= A[:, nodo] * valor

    # Imposición de Dirichlet
    for nodo, valor in nodos_dirichlet.items():

        A[nodo, :] = 0.0
        A[:, nodo] = 0.0

        A[nodo, nodo] = 1.0
        B[nodo] = valor

    return A, B


# ============================================================
# RESOLUCIÓN
# ============================================================

def resolver_fem_1d(L, Nelements, c, f, dirichlet_conditions, neumann_conditions):
    x, Nnodes = crear_malla(L, Nelements)

    A, B = ensamblar_sistema(x, Nelements, c, f)

    B = aplicar_neumann(B, neumann_conditions, Nnodes)

    A, B = aplicar_dirichlet(A, B, dirichlet_conditions, Nnodes)

    U = np.linalg.solve(A, B)

    U = U + 0.0

    return x, U, A, B


# ============================================================
# EJECUCIÓN
# ============================================================

x, U, A, B = resolver_fem_1d(
    L=L,
    Nelements=Nelements,
    c=c,
    f=f,
    dirichlet_conditions=dirichlet_conditions,
    neumann_conditions=neumann_conditions
)

np.set_printoptions(precision=10, suppress=True)

print("Matriz A:\n")

print(np.round(A, 10))

print("\nVector B:\n")

print(np.round(B, 10))

print("\nVector solución (u):\n")

print(np.round(U, 10))