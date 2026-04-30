import numpy as np
import meshio as ms
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from mpl_toolkits.mplot3d import Axes3D

# ============================================================
# LECTURA DE MALLA
# ============================================================

mesh = ms.read("C:/Users/santu/Desktop/Documentos ULPGC/Master/Asignaturas/Segundo Semestre/MSN/MSN/FEM 2D/pararayos.msh")

nodes = mesh.points[:, :2]
elements = mesh.cells_dict["triangle"]
lines = mesh.cells_dict["line"]
line_phys = mesh.cell_data_dict["gmsh:physical"]["line"]

Nnodes = len(nodes)
Nelements = len(elements)

print("Nodos:", Nnodes)
print("Triángulos:", Nelements)

# ============================================================
# PROBLEMA
# ============================================================

c = 1.0

def f(x, y):
    return 1.0

def g(x, y):
    return 1.0

# ============================================================
# MATRIZ GLOBAL
# ============================================================

A = np.zeros((Nnodes, Nnodes))
B = np.zeros(Nnodes)

# ============================================================
# CUADRATURA TRIANGULAR (3 puntos)
# ============================================================

gauss_pts = np.array([
    [1/6, 1/6],
    [2/3, 1/6],
    [1/6, 2/3]
])

weights = np.array([1/6, 1/6, 1/6])

# ============================================================
# ENSAMBLAJE
# ============================================================

grad_ref = np.array([
    [-1, -1],
    [ 1,  0],
    [ 0,  1]
])

for elem in elements:

    coords = nodes[elem]

    x1, y1 = coords[0]
    x2, y2 = coords[1]
    x3, y3 = coords[2]

    J = np.array([
        [x2 - x1, x3 - x1],
        [y2 - y1, y3 - y1]
    ])

    detJ = np.linalg.det(J)

    if abs(detJ) < 1e-14:
        raise ValueError("Elemento degenerado con detJ cercano a cero.")

    area_factor = abs(detJ)

    invJ = np.linalg.inv(J)

    # Transformación correcta de gradientes:
    # grad_phys = grad_ref @ invJ
    grad_phys = grad_ref @ invJ

    Ae = np.zeros((3, 3))
    Be = np.zeros(3)

    for k in range(3):

        xi, eta = gauss_pts[k]
        w = weights[k]

        N = np.array([
            1 - xi - eta,
            xi,
            eta
        ])

        x_phys = x1 + J[0, 0] * xi + J[0, 1] * eta
        y_phys = y1 + J[1, 0] * xi + J[1, 1] * eta

        for a in range(3):

            for b in range(3):

                # Término de rigidez: ∫ grad(Na) · grad(Nb) dΩ
                Ae[a, b] += w * np.dot(
                    grad_phys[a],
                    grad_phys[b]
                ) * area_factor

                # Término de masa/reacción: ∫ c Na Nb dΩ
                Ae[a, b] += w * c * N[a] * N[b] * area_factor

            # Segundo miembro: ∫ f Na dΩ
            Be[a] += w * f(x_phys, y_phys) * N[a] * area_factor

    # Ensamblaje global
    for a in range(3):
        A_global = elem[a]
        B[A_global] += Be[a]

        for b in range(3):
            B_global = elem[b]
            A[A_global, B_global] += Ae[a, b]

# ============================================================
# CONFIGURACIÓN DE CONDICIONES
# ============================================================

# Valores Dirichlet por grupo físico
dirichlet_values = {
    5: 0.0,
    6: 1.0
}

# Solo permitir 0 o 1 grupo Neumann
neumann_groups = []   # ejemplo: [7]

if len(neumann_groups) > 1:
    raise ValueError("Solo se permite una condición Neumann.")

# ============================================================
# DETECTAR NODOS DIRICHLET CON PRIORIDAD
# ============================================================

priority_groups = [6, 5]
priority_map = {group: i for i, group in enumerate(priority_groups)}

dirichlet_nodes = {}

for edge, phys in zip(lines, line_phys):

    if phys in dirichlet_values:

        value = dirichlet_values[phys]

        # Más robusto que priority_map[phys]
        priority = priority_map.get(phys, len(priority_groups))

        for node in edge:

            if node in dirichlet_nodes:

                _, old_priority = dirichlet_nodes[node]

                if priority < old_priority:
                    dirichlet_nodes[node] = (value, priority)

            else:
                dirichlet_nodes[node] = (value, priority)

dirichlet_nodes_set = set(dirichlet_nodes.keys())

# ============================================================
# APLICAR NEUMANN
# Dirichlet tiene prioridad
# ============================================================

if len(neumann_groups) == 1:

    neumann_group = neumann_groups[0]

    gauss_1D = np.array([
        0.5 - 1 / (2 * np.sqrt(3)),
        0.5 + 1 / (2 * np.sqrt(3))
    ])

    weights_1D = np.array([0.5, 0.5])

    for edge, phys in zip(lines, line_phys):

        if phys == neumann_group:

            n1, n2 = edge

            # Si los dos nodos son Dirichlet, no se añade nada
            if n1 in dirichlet_nodes_set and n2 in dirichlet_nodes_set:
                continue

            x1, y1 = nodes[n1]
            x2, y2 = nodes[n2]

            length = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)

            for xi, w in zip(gauss_1D, weights_1D):

                x_phys = x1 + xi * (x2 - x1)
                y_phys = y1 + xi * (y2 - y1)

                N = np.array([
                    1 - xi,
                    xi
                ])

                if n1 not in dirichlet_nodes_set:
                    B[n1] += w * g(x_phys, y_phys) * N[0] * length

                if n2 not in dirichlet_nodes_set:
                    B[n2] += w * g(x_phys, y_phys) * N[1] * length

# ============================================================
# APLICAR DIRICHLET CORRECTAMENTE
# ============================================================

fixed_nodes = list(dirichlet_nodes.keys())
fixed_values = np.array([dirichlet_nodes[node][0] for node in fixed_nodes])

# Corregir el segundo miembro antes de anular columnas
# Esto es necesario para condiciones Dirichlet no homogéneas, por ejemplo u = 1
for node, value in zip(fixed_nodes, fixed_values):
    B -= A[:, node] * value

# Imponer condiciones Dirichlet
for node, value in zip(fixed_nodes, fixed_values):

    A[node, :] = 0.0
    A[:, node] = 0.0
    A[node, node] = 1.0

    B[node] = value

# ============================================================
# RESOLVER
# ============================================================

U = np.linalg.solve(A, B)

print(f"Matriz A:\n{A}")
print(f"Tamaño de A: {A.shape}\n")

print(f"Vector B:\n{B}")
print(f"Tamaño de B: {B.shape}\n")

print(f"Solución U:\n{U}")
print(f"Tamaño de U: {U.shape}")

# ============================================================
# VISUALIZACIÓN
# ============================================================

triang = tri.Triangulation(nodes[:, 0], nodes[:, 1], elements)

# Malla
plt.figure()
plt.triplot(triang, color="black", linewidth=0.5)
plt.gca().set_aspect("equal")
plt.title("Malla FEM")
plt.xlabel("x")
plt.ylabel("y")
plt.show()

# Solución 2D
plt.figure()
plt.tricontourf(triang, U, 50)
plt.colorbar(label="u(x,y)")
plt.gca().set_aspect("equal")
plt.title("Solución FEM")
plt.xlabel("x")
plt.ylabel("y")
plt.show()

# Solución 3D
fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")
ax.plot_trisurf(
    nodes[:, 0],
    nodes[:, 1],
    U,
    triangles=elements,
    cmap="viridis"
)
ax.set_title("Superficie FEM")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("u")
plt.show()