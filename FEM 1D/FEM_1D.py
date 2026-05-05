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

# def u_real(x):
#     return ...

u_real = None


# ============================================================
# MALLA
# ============================================================

def crear_malla(L, Nelements):
    Nnodes = Nelements + 1
    x = np.linspace(0.0, L, Nnodes)
    return x, Nnodes


# ============================================================
# CUADRATURA Y FUNCIONES DE FORMA EN [0, 1]
# ============================================================

def gauss_2p():
    xi = np.array([
        (1.0 - 1.0 / np.sqrt(3.0)) / 2.0,
        (1.0 + 1.0 / np.sqrt(3.0)) / 2.0
    ])

    w = np.array([
        0.5,
        0.5
    ])

    return xi, w


def funciones_forma_lineales(xi):
    N = np.array([
        1.0 - xi,
        xi
    ])

    dN_dxi = np.array([
        -1.0,
         1.0
    ])

    return N, dN_dxi


# ============================================================
# MATRICES LOCALES
# ============================================================

def calcular_matriz_local(x1, x2, c, f):
    h = x2 - x1
    J = h

    Ke = (1.0 / h) * np.array([
        [ 1.0, -1.0],
        [-1.0,  1.0]
    ])

    Me = np.zeros((2, 2))

    Be = np.zeros(2)

    xi_gauss, pesos = gauss_2p()

    for xi, w in zip(xi_gauss, pesos):
        N, _ = funciones_forma_lineales(xi)

        # Transformación al elemento físico:
        x_phys = x1 + h * xi

        for a in range(2):

            for b in range(2):

                Me[a, b] += (
                    w
                    * N[a]
                    * N[b]
                    * J
                )

            Be[a] += (
                w
                * f(x_phys)
                * N[a]
                * J
            )

    Ae = Ke + c * Me

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

    for nodo, valor in nodos_dirichlet.items():
        B -= A[:, nodo] * valor

    for nodo, valor in nodos_dirichlet.items():

        A[nodo, :] = 0.0
        A[:, nodo] = 0.0

        A[nodo, nodo] = 1.0
        B[nodo] = valor

    return A, B


# ============================================================
# RESOLUCIÓN
# ============================================================

def resolver_fem_1d(
    L,
    Nelements,
    c,
    f,
    dirichlet_conditions,
    neumann_conditions
):
    x, Nnodes = crear_malla(L, Nelements)

    A, B = ensamblar_sistema(x, Nelements, c, f)

    B = aplicar_neumann(
        B,
        neumann_conditions,
        Nnodes
    )

    A, B = aplicar_dirichlet(
        A,
        B,
        dirichlet_conditions,
        Nnodes
    )

    U = np.linalg.solve(A, B)

    return x, U, A, B


# ============================================================
# CÁLCULO DEL ERROR
# ============================================================

def calcular_error(x, U, u_exacta):
    if u_exacta is None:
        return None

    U_exacta = np.array([
        u_exacta(xi) for xi in x
    ])

    error = U_exacta - U

    error_max = np.max(np.abs(error))

    return U_exacta, error, error_max


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

resultado_error = calcular_error(x, U, u_real)

np.set_printoptions(precision=10, suppress=True)

print("\nMatriz A:")
print(np.round(A, 10))
print(f"Tamaño de A: {A.shape}\n")

print("\nVector B:")
print(np.round(B, 10))
print(f"Tamaño de B: {B.shape}\n")

print("\nVector solución U:")
print(np.round(U, 10))
print(f"Tamaño de U: {U.shape}\n")

if resultado_error is not None:
    _, error, error_max = resultado_error

    print("\nError U_exacta - U_FEM:")
    print(np.round(error, 10))

    print("\nError máximo:")
    print(error_max)
else:
    print("\nNo se ha definido una solución exacta, por lo que no se puede calcular el error.")
