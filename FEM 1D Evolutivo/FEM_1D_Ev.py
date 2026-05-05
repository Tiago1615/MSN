import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


# ============================================================
# CONFIGURACIÓN DEL PROBLEMA
# ============================================================

L = 10.0
Nelements = 100

alpha = 1.0
dt = 0.01
Tfinal = 1.0

guardar_txt = True
guardar_graficas = True
guardar_animacion = True


def f(x, t):
    return np.cos(2.0 * np.pi * x)


def u_inicial(x):
    return 0.0


dirichlet_conditions = {
    "left": 0.0,
    "right": 1.0
}

neumann_conditions = {}

u_exacta = None

# def u_exacta(x, t):
#     ...


# ============================================================
# MALLA
# ============================================================

def crear_malla(L, Nelements):
    Nnodes = Nelements + 1
    x = np.linspace(0.0, L, Nnodes)
    return x


# ============================================================
# CUADRATURA Y FUNCIONES DE FORMA EN [0, 1]
# ============================================================

def cuadratura_gauss_2p():
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
# UTILIDADES
# ============================================================

def nodo_lado(lado, Nnodes):
    if lado == "left":
        return 0

    if lado == "right":
        return Nnodes - 1

    raise ValueError(f"Lado no reconocido: {lado}")


def evaluar_valor(valor, x, t):
    if callable(valor):
        return valor(x, t)

    return float(valor)


def aplicar_dirichlet_a_vector(U, x, dirichlet_conditions, t):
    Nnodes = len(x)

    for lado, valor in dirichlet_conditions.items():
        nodo = nodo_lado(lado, Nnodes)

        U[nodo] = evaluar_valor(
            valor,
            x[nodo],
            t
        )

    return U


# ============================================================
# MATRICES LOCALES
# ============================================================

def matrices_locales_1d(x1, x2):
    h = x2 - x1
    J = h

    Ke = (1.0 / h) * np.array([
        [ 1.0, -1.0],
        [-1.0,  1.0]
    ])

    Me = np.zeros((2, 2))

    xi_gauss, pesos = cuadratura_gauss_2p()

    for xi, w in zip(xi_gauss, pesos):

        N, _ = funciones_forma_lineales(xi)

        for a in range(2):

            for b in range(2):

                Me[a, b] += (
                    w
                    * N[a]
                    * N[b]
                    * J
                )

    return Ke, Me


def vector_cargas_local_1d(x1, x2, f, t):
    h = x2 - x1
    J = h

    Fe = np.zeros(2)

    xi_gauss, pesos = cuadratura_gauss_2p()

    for xi, w in zip(xi_gauss, pesos):

        N, _ = funciones_forma_lineales(xi)

        # Transformación:
        x_phys = x1 + h * xi

        for a in range(2):

            Fe[a] += (
                w * f(x_phys, t) * N[a] * J
            )

    return Fe


# ============================================================
# ENSAMBLAJE GLOBAL
# ============================================================

def ensamblar_matrices(x):
    Nnodes = len(x)
    Nelements = Nnodes - 1

    K = np.zeros((Nnodes, Nnodes))
    M = np.zeros((Nnodes, Nnodes))

    for e in range(Nelements):

        x1 = x[e]
        x2 = x[e + 1]

        Ke, Me = matrices_locales_1d(x1, x2)

        nodos = [e, e + 1]

        for a in range(2):

            I = nodos[a]

            for b in range(2):

                Jg = nodos[b]

                K[I, Jg] += Ke[a, b]
                M[I, Jg] += Me[a, b]

    return K, M


def ensamblar_vector_cargas(x, f, t):
    Nnodes = len(x)
    Nelements = Nnodes - 1

    F = np.zeros(Nnodes)

    for e in range(Nelements):

        x1 = x[e]
        x2 = x[e + 1]

        Fe = vector_cargas_local_1d(x1, x2, f, t)

        nodos = [e, e + 1]

        for a in range(2):

            I = nodos[a]

            F[I] += Fe[a]

    return F


# ============================================================
# CONDICIONES DE CONTORNO
# ============================================================

def aplicar_neumann_rhs(b, x, neumann_conditions, alpha, dt, t):
    Nnodes = len(x)

    for lado, valor in neumann_conditions.items():

        nodo = nodo_lado(lado, Nnodes)

        g = evaluar_valor(
            valor,
            x[nodo],
            t
        )

        if lado == "left":
            b[nodo] -= dt * alpha * g

        elif lado == "right":
            b[nodo] += dt * alpha * g

    return b


def aplicar_dirichlet_sistema(A, b, x, dirichlet_conditions, t):
    Nnodes = len(x)

    dirichlet_nodes = {}

    for lado, valor in dirichlet_conditions.items():

        nodo = nodo_lado(lado, Nnodes)

        valor_nodo = evaluar_valor(
            valor,
            x[nodo],
            t
        )

        dirichlet_nodes[nodo] = valor_nodo

    for nodo, valor in dirichlet_nodes.items():
        b -= A[:, nodo] * valor

    # Imponer Dirichlet
    for nodo, valor in dirichlet_nodes.items():

        A[nodo, :] = 0.0
        A[:, nodo] = 0.0

        A[nodo, nodo] = 1.0
        b[nodo] = valor

    return A, b


# ============================================================
# RESOLUCIÓN TEMPORAL
# ============================================================

def resolver_fem_1d_evolutivo(
    L,
    Nelements,
    alpha,
    dt,
    Tfinal,
    f,
    u_inicial,
    dirichlet_conditions,
    neumann_conditions
):
    x = crear_malla(L, Nelements)

    Nt = int(Tfinal / dt)

    K, M = ensamblar_matrices(x)

    # Condición inicial
    U = np.array(
        [u_inicial(xi) for xi in x],
        dtype=float
    )

    U = aplicar_dirichlet_a_vector(
        U,
        x,
        dirichlet_conditions,
        t=0.0
    )

    U_history = [U.copy()]
    time_history = [0.0]

    A_base = M + dt * alpha * K

    for n in range(Nt):
        t_next = (n + 1) * dt

        F = ensamblar_vector_cargas(
            x,
            f,
            t_next
        )

        b = M @ U + dt * F

        b = aplicar_neumann_rhs(
            b,
            x,
            neumann_conditions,
            alpha,
            dt,
            t_next
        )

        A_step = A_base.copy()

        A_step, b = aplicar_dirichlet_sistema(
            A_step,
            b,
            x,
            dirichlet_conditions,
            t_next
        )

        U = np.linalg.solve(A_step, b)

        U_history.append(U.copy())
        time_history.append(t_next)

    U_array = np.array(U_history)
    time_array = np.array(time_history)

    return x, time_array, U_array, K, M


# ============================================================
# GUARDAR RESULTADOS
# ============================================================

def guardar_soluciones_txt(x, time_array, U_array, carpeta):
    os.makedirs(carpeta, exist_ok=True)

    for i, (t, U) in enumerate(zip(time_array, U_array)):

        nombre_archivo = os.path.join(
            carpeta,
            f"solucion_t_{i:04d}.txt"
        )

        datos = np.column_stack((x, U))

        np.savetxt(
            nombre_archivo,
            datos,
            header=f"t = {t:.6f}\nx    u(x,t)",
            comments=""
        )

    print(f"Soluciones temporales guardadas en la carpeta: {carpeta}")


# ============================================================
# VISUALIZACIÓN
# ============================================================

def graficar_curvas(x, time_array, U_array, archivo=None):
    Nt = len(time_array) - 1

    indices = [
        0,
        Nt // 4,
        Nt // 2,
        3 * Nt // 4,
        Nt
    ]

    plt.figure(figsize=(8, 5))

    for idx in indices:
        plt.plot(
            x,
            U_array[idx],
            label=f"t={time_array[idx]:.2f}"
        )

    plt.xlabel("x")
    plt.ylabel("u(x,t)")
    plt.title("Evolución temporal FEM 1D")
    plt.legend()
    plt.grid(True)

    if archivo is not None:
        plt.savefig(archivo)

    plt.show()


def graficar_mapa_calor(x, time_array, U_array, archivo=None):
    plt.figure(figsize=(8, 5))

    plt.imshow(
        U_array,
        aspect="auto",
        extent=[
            x[0],
            x[-1],
            time_array[0],
            time_array[-1]
        ],
        origin="lower"
    )

    plt.colorbar(label="u(x,t)")
    plt.xlabel("x")
    plt.ylabel("t")
    plt.title("Mapa de calor")

    if archivo is not None:
        plt.savefig(archivo)

    plt.show()


def guardar_animacion(x, time_array, U_array, archivo="animacion.gif"):
    fig, ax = plt.subplots(figsize=(8, 5))

    line, = ax.plot(
        x,
        U_array[0],
        lw=2
    )

    ax.set_xlim(x[0], x[-1])

    ymin = np.min(U_array)
    ymax = np.max(U_array)

    margen = 0.05 * (ymax - ymin) if ymax > ymin else 0.1

    ax.set_ylim(
        ymin - margen,
        ymax + margen
    )

    ax.set_xlabel("x")
    ax.set_ylabel("u(x,t)")
    ax.grid(True)

    def update(frame):
        line.set_ydata(U_array[frame])
        ax.set_title(f"t = {time_array[frame]:.3f}")
        return line,

    ani = FuncAnimation(
        fig,
        update,
        frames=len(time_array),
        interval=50,
        blit=True
    )

    ani.save(
        archivo,
        writer="pillow",
        fps=20
    )

    plt.close(fig)

    print(f"Animación guardada como {archivo}")


# ============================================================
# CÁLCULO DEL ERROR
# ============================================================

def calcular_error(x, time_array, U_array, u_exacta):
    if u_exacta is None:
        return None

    errores_max = []

    for n, t in enumerate(time_array):

        U = U_array[n]

        U_exacta = np.array([
            u_exacta(xi, t) for xi in x
        ])

        error = U_exacta - U

        errores_max.append(
            np.max(np.abs(error))
        )

    return U_exacta, np.array(error), np.array(errores_max)


# ============================================================
# EJECUCIÓN
# ============================================================

x, time_array, U_array, K, M = resolver_fem_1d_evolutivo(
    L=L,
    Nelements=Nelements,
    alpha=alpha,
    dt=dt,
    Tfinal=Tfinal,
    f=f,
    u_inicial=u_inicial,
    dirichlet_conditions=dirichlet_conditions,
    neumann_conditions=neumann_conditions
)

resultado_error = calcular_error(
    x,
    time_array,
    U_array,
    u_exacta
)

np.set_printoptions(precision=10, suppress=True)

print("\nMatriz K:")
print(np.round(K, 10))
print(f"Tamaño de K: {K.shape}\n")

print("\nMatriz M:")
print(np.round(M, 10))
print(f"Tamaño de M: {M.shape}\n")

print(f"\nValor del vector solución en t= {Tfinal:.2f}:")
print(U_array[-1])
print(f"Tamaño de U_array: {U_array.shape}\n")

if resultado_error is not None:
    _, error, errores_max = resultado_error

    print("\nError U_exacta - U_FEM en el tiempo final:")
    print(np.round(error, 10))

    print("\nError máximo en cada tiempo:")
    print(errores_max)
else:
    print("\nNo se ha definido una función de solución exacta, por lo que no se calcula el error.\n")

if guardar_txt:
    guardar_soluciones_txt(
        x,
        time_array,
        U_array,
        carpeta="soluciones_temporales"
    )

if guardar_graficas:
    graficar_curvas(
        x,
        time_array,
        U_array,
        archivo="curvas.png"
    )

    graficar_mapa_calor(
        x,
        time_array,
        U_array,
        archivo="heatmap.png"
    )

if guardar_animacion:
    guardar_animacion(
        x,
        time_array,
        U_array,
        archivo="animacion.gif"
    )
