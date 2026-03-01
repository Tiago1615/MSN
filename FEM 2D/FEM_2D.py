import numpy as np
import meshio as ms
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from mpl_toolkits.mplot3d import Axes3D

# ============================================================
# 1) LECTURA DE MALLA
# ============================================================

mesh = ms.read("C:/Users/santu/Desktop/Documentos ULPGC/Master/Asignaturas/Segundo Semestre/MSN/MSN/FEM 2D/cuadrado.msh")

nodes = mesh.points[:, :2]
elements = mesh.cells_dict["triangle"]
lines = mesh.cells_dict["line"]
line_phys = mesh.cell_data_dict["gmsh:physical"]["line"]

Nnodes = len(nodes)
Nelements = len(elements)

print("Nodos:", Nnodes)
print("Triángulos:", Nelements)

# ============================================================
# 2) PROBLEMA
# ============================================================

c = 1.0

def f(x, y):
    return 0.0

def g(x, y):
    return 1.0   # valor Neumann si se usa

# ============================================================
# 3) MATRIZ GLOBAL
# ============================================================

A = np.zeros((Nnodes, Nnodes))
B = np.zeros(Nnodes)

# ============================================================
# 4) CUADRATURA TRIANGULAR (3 puntos)
# ============================================================

gauss_pts = np.array([
    [1/6, 1/6],
    [2/3, 1/6],
    [1/6, 2/3]
])

weights = np.array([1/6, 1/6, 1/6])

# ============================================================
# 5) ENSAMBLAJE
# ============================================================

for elem in elements:

    coords = nodes[elem]

    x1, y1 = coords[0]
    x2, y2 = coords[1]
    x3, y3 = coords[2]

    J = np.array([
        [x2-x1, x3-x1],
        [y2-y1, y3-y1]
    ])

    detJ = np.linalg.det(J)
    invJ = np.linalg.inv(J)

    grad_ref = np.array([
        [-1, -1],
        [ 1,  0],
        [ 0,  1]
    ])

    grad_phys = grad_ref @ invJ.T

    Ae = np.zeros((3,3))
    Be = np.zeros(3)

    for k in range(3):

        xi, eta = gauss_pts[k]
        w = weights[k]

        N = np.array([1-xi-eta, xi, eta])

        x_phys = x1 + J[0,0]*xi + J[0,1]*eta
        y_phys = y1 + J[1,0]*xi + J[1,1]*eta

        for a in range(3):
            for b in range(3):

                # Rigidez
                Ae[a,b] += w * (
                    np.dot(grad_phys[a], grad_phys[b]) * detJ
                )

                # Masa
                Ae[a,b] += w * (
                    c * N[a] * N[b] * detJ
                )

            Be[a] += w * (
                f(x_phys, y_phys) * N[a] * detJ
            )

    for a in range(3):
        for b in range(3):
            A[elem[a], elem[b]] += Ae[a,b]
        B[elem[a]] += Be[a]

# ============================================================
# 6) CONFIGURACIÓN DE CONDICIONES
# ============================================================

# Valores Dirichlet por grupo físico
dirichlet_values = {
    5: 0.0,
    6: 1.0
}

# Solo permitir 0 o 1 grupo Neumann
neumann_groups = []   # ejemplo: [6]

if len(neumann_groups) > 1:
    raise ValueError("Solo se permite una condición Neumann.")

# ============================================================
# 7) DETECTAR NODOS DIRICHLET CON PRIORIDAD
# ============================================================

priority_groups = [6, 5]
priority_map = {group: i for i, group in enumerate(priority_groups)}

dirichlet_nodes = {}

for edge, phys in zip(lines, line_phys):

    if phys in dirichlet_values:

        value = dirichlet_values[phys]
        priority = priority_map[phys]

        for node in edge:

            if node in dirichlet_nodes:
                _, old_priority = dirichlet_nodes[node]
                if priority < old_priority:
                    dirichlet_nodes[node] = (value, priority)
            else:
                dirichlet_nodes[node] = (value, priority)

dirichlet_nodes_set = set(dirichlet_nodes.keys())

# ============================================================
# 8) APLICAR NEUMANN (Dirichlet tiene prioridad)
# ============================================================

if len(neumann_groups) == 1:

    neumann_group = neumann_groups[0]

    gauss_1D = np.array([
        0.5 - 1/(2*np.sqrt(3)),
        0.5 + 1/(2*np.sqrt(3))
    ])

    weights_1D = np.array([0.5, 0.5])

    for edge, phys in zip(lines, line_phys):

        if phys == neumann_group:

            n1, n2 = edge

            if n1 in dirichlet_nodes_set and n2 in dirichlet_nodes_set:
                continue

            x1, y1 = nodes[n1]
            x2, y2 = nodes[n2]

            length = np.sqrt((x2-x1)**2 + (y2-y1)**2)

            for xi, w in zip(gauss_1D, weights_1D):

                x_phys = x1 + xi*(x2-x1)
                y_phys = y1 + xi*(y2-y1)

                N = np.array([1-xi, xi])

                if n1 not in dirichlet_nodes_set:
                    B[n1] += w * g(x_phys,y_phys) * N[0] * length

                if n2 not in dirichlet_nodes_set:
                    B[n2] += w * g(x_phys,y_phys) * N[1] * length

# ============================================================
# 9) APLICAR DIRICHLET
# ============================================================

for node, (value, _) in dirichlet_nodes.items():

    A[node,:] = 0.0
    A[:,node] = 0.0
    A[node,node] = 1.0
    B[node] = value

# ============================================================
# 10) RESOLVER
# ============================================================

U = np.linalg.solve(A, B)

print("Sistema resuelto correctamente.")

# ============================================================
# 11) VISUALIZACIÓN
# ============================================================

triang = tri.Triangulation(nodes[:,0], nodes[:,1], elements)

# Malla
plt.figure()
plt.triplot(triang, color="black", linewidth=0.5)
plt.gca().set_aspect("equal")
plt.title("Malla FEM")
plt.show()

# Solución 2D
plt.figure()
plt.tricontourf(triang, U, 50)
plt.colorbar(label="u(x,y)")
plt.gca().set_aspect("equal")
plt.title("Solución FEM")
plt.show()

# Solución 3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_trisurf(nodes[:,0], nodes[:,1], U, triangles=elements, cmap='viridis')
ax.set_title("Superficie FEM")
plt.show()