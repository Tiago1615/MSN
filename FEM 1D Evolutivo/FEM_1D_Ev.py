import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# ============================================================
# PARÁMETROS DEL PROBLEMA
# ============================================================

L = 10.0
Nelements = 100
Nnodes = Nelements + 1

alpha = 1.0
dt = 0.01
Tfinal = 1.0
Nt = int(Tfinal / dt)

def f(x):
    return np.cos(2.0 * np.pi * x)

# ============================================================
# MALLA
# ============================================================

x = np.linspace(0.0, L, Nnodes)

# ============================================================
# MATRICES GLOBALES
# ============================================================

K = np.zeros((Nnodes, Nnodes))
M = np.zeros((Nnodes, Nnodes))
F = np.zeros(Nnodes)

# ============================================================
# CUADRATURA DE GAUSS EN [-1, 1]
# ============================================================

xi_gauss = np.array([
    -1.0 / np.sqrt(3.0),
     1.0 / np.sqrt(3.0)
])

weights = np.array([1.0, 1.0])

# ============================================================
# ENSAMBLAJE FEM 1D
# ============================================================

for e in range(Nelements):

    x1 = x[e]
    x2 = x[e + 1]
    h = x2 - x1

    Ke = np.zeros((2, 2))
    Me = np.zeros((2, 2))
    Fe = np.zeros(2)

    for k in range(2):

        xi = xi_gauss[k]
        w = weights[k]

        # Funciones de forma en [-1, 1]
        Nhat = np.array([
            (1.0 - xi) / 2.0,
            (1.0 + xi) / 2.0
        ])

        dNhat_dxi = np.array([
            -0.5,
             0.5
        ])

        # Mapeo a coordenada física
        x_phys = (x1 + x2) / 2.0 + (h / 2.0) * xi

        # Jacobiano
        J = h / 2.0

        # Derivadas respecto a x
        dNdx = dNhat_dxi * (2.0 / h)

        for i in range(2):

            for j in range(2):

                # Matriz de rigidez local
                Ke[i, j] += w * dNdx[i] * dNdx[j] * J

                # Matriz de masa local
                Me[i, j] += w * Nhat[i] * Nhat[j] * J

            # Vector de cargas local
            Fe[i] += w * f(x_phys) * Nhat[i] * J

    # Ensamblaje global
    nodes = [e, e + 1]

    for i in range(2):

        I = nodes[i]
        F[I] += Fe[i]

        for j in range(2):

            Jg = nodes[j]
            K[I, Jg] += Ke[i, j]
            M[I, Jg] += Me[i, j]

# ============================================================
# CONDICIONES DE DIRICHLET
# ============================================================

def apply_dirichlet(A, b, dirichlet_values):
    """
    Aplica condiciones de Dirichlet de forma simétrica.

    dirichlet_values: diccionario {nodo: valor}
    """

    # Primero corregimos el segundo miembro
    for node, value in dirichlet_values.items():
        b -= A[:, node] * value

    # Luego anulamos filas y columnas
    for node, value in dirichlet_values.items():

        A[node, :] = 0.0
        A[:, node] = 0.0

        A[node, node] = 1.0
        b[node] = value

# ============================================================
# CONDICIÓN INICIAL
# ============================================================

U = np.zeros(Nnodes)

# Aplicamos Dirichlet desde el inicio
U[0] = 0.0
U[-1] = 1.0

dirichlet_values = {
    0: 0.0,
    Nnodes - 1: 1.0
}

U_history = [U.copy()]
time_history = [0.0]

# ============================================================
# BUCLE TEMPORAL
# ============================================================

# Esquema de Euler implícito:
#
# (M + dt * alpha * K) U^{n+1} = M U^n + dt F

A_base = M + dt * alpha * K

for n in range(Nt):

    b = M @ U + dt * F

    A_step = A_base.copy()

    apply_dirichlet(A_step, b, dirichlet_values)

    U = np.linalg.solve(A_step, b)

    U_history.append(U.copy())
    time_history.append((n + 1) * dt)

U_array = np.array(U_history)
time_array = np.array(time_history)

print(f"Valor del vector solución en t={Tfinal:.2f}:")
print(U)

print(f"\nTamaño de U_array: {U_array.shape}")

# ============================================================
# GUARDAR SOLUCIONES TEMPORALES EN TXT
# Igual que en el segundo código
# ============================================================

import os

carpeta = "soluciones_temporales"
os.makedirs(carpeta, exist_ok=True)

for i, (sol, t) in enumerate(zip(U_array, time_array)):

    nombre_archivo = os.path.join(carpeta, f"solucion_t_{i:04d}.txt")

    datos = np.column_stack((x, sol))

    np.savetxt(
        nombre_archivo,
        datos,
        header=f"t = {t:.6f}\nx    u(x,t)",
        comments=""
    )

print(f"\nSoluciones temporales guardadas en la carpeta: {carpeta}")

# ============================================================
# GRÁFICAS
# ============================================================

# ------------------------------------------------------------
# Curvas en distintos tiempos
# ------------------------------------------------------------

plt.figure(figsize=(8, 5))

times_to_plot = [
    0,
    Nt // 4,
    Nt // 2,
    3 * Nt // 4,
    Nt
]

for t_idx in times_to_plot:
    plt.plot(
        x,
        U_array[t_idx],
        label=f"t={time_array[t_idx]:.2f}"
    )

plt.xlabel("x")
plt.ylabel("u(x,t)")
plt.title("Evolución temporal FEM 1D")
plt.legend()
plt.grid(True)

plt.savefig("curvas.png")
plt.show()

# ------------------------------------------------------------
# Mapa de calor
# ------------------------------------------------------------

plt.figure(figsize=(8, 5))

plt.imshow(
    U_array,
    aspect="auto",
    extent=[x[0], x[-1], time_array[0], time_array[-1]],
    origin="lower"
)

plt.colorbar(label="u(x,t)")
plt.xlabel("x")
plt.ylabel("t")
plt.title("Mapa de calor")

plt.savefig("heatmap.png")
plt.show()

# ------------------------------------------------------------
# Animación
# ------------------------------------------------------------

fig, ax = plt.subplots(figsize=(8, 5))

line, = ax.plot(x, U_array[0], lw=2)

ax.set_xlim(x[0], x[-1])
ax.set_ylim(np.min(U_array), np.max(U_array))
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
    frames=len(U_array),
    interval=50,
    blit=True
)

ani.save("animacion.gif", writer="pillow", fps=20)

plt.close(fig)

print("Animación guardada como animacion.gif")