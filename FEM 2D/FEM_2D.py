import numpy as np
import meshio as ms
import matplotlib.pyplot as plt
import matplotlib.tri as tri


# ============================================================
# CONFIGURACIÓN DEL PROBLEMA
# ============================================================

mesh_file = "C:/Users/santu/Desktop/Documentos ULPGC/Master/Asignaturas/Segundo Semestre/MSN/MSN/FEM 2D/pararayos.msh"

c = 1.0


def f(x, y):
    return 1.0


def g(x, y):
    return 1.0

dirichlet_conditions = {
    5: 0.0,
    6: 1.0
}

neumann_conditions = {}

priority_groups = [6, 5]


# ============================================================
# LECTURA DE MALLA
# ============================================================

def leer_malla(archivo):
    mesh = ms.read(archivo)

    nodes = mesh.points[:, :2]
    elements = mesh.cells_dict["triangle"]
    lines = mesh.cells_dict["line"]
    line_phys = mesh.cell_data_dict["gmsh:physical"]["line"]

    return nodes, elements, lines, line_phys


# ============================================================
# CUADRATURA
# ============================================================

def cuadratura_triangulo_3p():
    puntos = np.array([
        [1.0 / 6.0, 1.0 / 6.0],
        [2.0 / 3.0, 1.0 / 6.0],
        [1.0 / 6.0, 2.0 / 3.0]
    ])

    pesos = np.array([
        1.0 / 6.0,
        1.0 / 6.0,
        1.0 / 6.0
    ])

    return puntos, pesos


def cuadratura_linea_2p():
    puntos = np.array([
        (1.0 - 1.0 / np.sqrt(3.0)) / 2.0,
        (1.0 + 1.0 / np.sqrt(3.0)) / 2.0
    ])

    pesos = np.array([
        0.5,
        0.5
    ])

    return puntos, pesos


# ============================================================
# FUNCIONES DE FORMA
# ============================================================

def funciones_forma_triangulo(xi, eta):
    return np.array([
        1.0 - xi - eta,
        xi,
        eta
    ])


def gradientes_referencia():
    return np.array([
        [-1.0, -1.0],
        [ 1.0,  0.0],
        [ 0.0,  1.0]
    ])


def funciones_forma_linea(xi):
    return np.array([
        1.0 - xi,
        xi
    ])


# ============================================================
# UTILIDADES
# ============================================================

def evaluar_valor(valor, x, y):
    if callable(valor):
        return valor(x, y)

    return float(valor)


def matriz_jacobiana(coords):
    x1, y1 = coords[0]
    x2, y2 = coords[1]
    x3, y3 = coords[2]

    J = np.array([
        [x2 - x1, x3 - x1],
        [y2 - y1, y3 - y1]
    ])

    return J


def mapear_a_fisico(coords, xi, eta):
    x1, y1 = coords[0]

    J = matriz_jacobiana(coords)

    x = x1 + J[0, 0] * xi + J[0, 1] * eta
    y = y1 + J[1, 0] * xi + J[1, 1] * eta

    return x, y


def gradientes_fisicos(coords):
    J = matriz_jacobiana(coords)

    invJT = np.linalg.inv(J).T

    grad_ref = gradientes_referencia()

    grad_phys = np.zeros((3, 2))

    for a in range(3):
        grad_phys[a, :] = invJT @ grad_ref[a, :]

    return grad_phys


# ============================================================
# MATRICES LOCALES
# ============================================================

def calcular_elemento_triangular(coords, c, f):
    J = matriz_jacobiana(coords)

    detJ = np.linalg.det(J)

    if abs(detJ) < 1e-14:
        raise ValueError("Elemento degenerado con detJ cercano a cero.")

    area_factor = abs(detJ)
    area_K = 0.5 * abs(detJ)

    grad_phys = gradientes_fisicos(coords)

    Ke = np.zeros((3, 3))

    for a in range(3):
        for b in range(3):
            Ke[a, b] = (
                np.dot(grad_phys[a], grad_phys[b])
                * area_K
            )

    Me = np.zeros((3, 3))
    Be = np.zeros(3)

    puntos, pesos = cuadratura_triangulo_3p()

    for (xi, eta), w in zip(puntos, pesos):

        N = funciones_forma_triangulo(xi, eta)

        x_phys, y_phys = mapear_a_fisico(coords, xi, eta)

        for a in range(3):

            for b in range(3):

                Me[a, b] += (
                    w
                    * N[a]
                    * N[b]
                    * area_factor
                )

            Be[a] += (
                w
                * f(x_phys, y_phys)
                * N[a]
                * area_factor
            )

    Ae = Ke + c * Me
    return Ae, Be


# ============================================================
# ENSAMBLAJE GLOBAL
# ============================================================

def ensamblar_sistema(nodes, elements, c, f):
    Nnodes = len(nodes)

    A = np.zeros((Nnodes, Nnodes))
    B = np.zeros(Nnodes)

    for elem in elements:

        coords = nodes[elem]

        Ae, Be = calcular_elemento_triangular(coords, c, f)

        for a in range(3):

            I = elem[a]
            B[I] += Be[a]

            for b in range(3):

                Jg = elem[b]
                A[I, Jg] += Ae[a, b]

    return A, B


# ============================================================
# CONDICIONES DE DIRICHLET
# ============================================================

def detectar_nodos_dirichlet(
    nodes,
    lines,
    line_phys,
    dirichlet_conditions,
    priority_groups
):
    priority_map = {
        group: i for i, group in enumerate(priority_groups)
    }

    dirichlet_nodes = {}

    for edge, phys in zip(lines, line_phys):

        if phys not in dirichlet_conditions:
            continue

        condicion = dirichlet_conditions[phys]

        priority = priority_map.get(
            phys,
            len(priority_groups)
        )

        for node in edge:

            x_node, y_node = nodes[node]

            value = evaluar_valor(
                condicion,
                x_node,
                y_node
            )

            if node in dirichlet_nodes:

                _, old_priority = dirichlet_nodes[node]

                if priority < old_priority:
                    dirichlet_nodes[node] = (value, priority)

            else:
                dirichlet_nodes[node] = (value, priority)

    resultado = {
        node: value_priority[0]
        for node, value_priority in dirichlet_nodes.items()
    }

    return resultado


def aplicar_dirichlet(A, B, dirichlet_nodes):
    for node, value in dirichlet_nodes.items():
        B -= A[:, node] * value

    for node, value in dirichlet_nodes.items():

        A[node, :] = 0.0
        A[:, node] = 0.0

        A[node, node] = 1.0
        B[node] = value

    return A, B


# ============================================================
# CONDICIONES DE NEUMANN
# ============================================================

def aplicar_neumann(
    B,
    nodes,
    lines,
    line_phys,
    neumann_conditions,
    dirichlet_nodes
):
    dirichlet_set = set(dirichlet_nodes.keys())

    puntos, pesos = cuadratura_linea_2p()

    for edge, phys in zip(lines, line_phys):

        if phys not in neumann_conditions:
            continue

        condicion = neumann_conditions[phys]

        n1, n2 = edge

        x1, y1 = nodes[n1]
        x2, y2 = nodes[n2]

        length = np.sqrt(
            (x2 - x1)**2
            + (y2 - y1)**2
        )

        for xi, w in zip(puntos, pesos):

            x_phys = x1 + xi * (x2 - x1)
            y_phys = y1 + xi * (y2 - y1)

            N = funciones_forma_linea(xi)

            g_value = evaluar_valor(
                condicion,
                x_phys,
                y_phys
            )

            if n1 not in dirichlet_set:
                B[n1] += w * g_value * N[0] * length

            if n2 not in dirichlet_set:
                B[n2] += w * g_value * N[1] * length

    return B


# ============================================================
# RESOLUCIÓN FEM 2D
# ============================================================

def resolver_fem_2d(
    mesh_file,
    c,
    f,
    dirichlet_conditions,
    neumann_conditions,
    priority_groups
):
    nodes, elements, lines, line_phys = leer_malla(mesh_file)

    print("Número de Nodos:", len(nodes))
    print("Número de Triángulos:", len(elements))

    A, B = ensamblar_sistema(
        nodes,
        elements,
        c,
        f
    )

    dirichlet_nodes = detectar_nodos_dirichlet(
        nodes,
        lines,
        line_phys,
        dirichlet_conditions,
        priority_groups
    )

    B = aplicar_neumann(
        B,
        nodes,
        lines,
        line_phys,
        neumann_conditions,
        dirichlet_nodes
    )

    A, B = aplicar_dirichlet(
        A,
        B,
        dirichlet_nodes
    )

    U = np.linalg.solve(A, B)

    return nodes, elements, A, B, U


# ============================================================
# VISUALIZACIÓN
# ============================================================

def visualizar_resultado(nodes, elements, U):
    triang = tri.Triangulation(
        nodes[:, 0],
        nodes[:, 1],
        elements
    )

    # Malla
    plt.figure()
    plt.triplot(
        triang,
        color="black",
        linewidth=0.5
    )
    plt.gca().set_aspect("equal")
    plt.title("Malla FEM")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.show()

    # Solución 2D
    plt.figure()
    plt.tricontourf(
        triang,
        U,
        50
    )
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


# ============================================================
# EJECUCIÓN
# ============================================================

nodes, elements, A, B, U = resolver_fem_2d(
    mesh_file=mesh_file,
    c=c,
    f=f,
    dirichlet_conditions=dirichlet_conditions,
    neumann_conditions=neumann_conditions,
    priority_groups=priority_groups
)

np.set_printoptions(precision=10, suppress=True)

print(f"\nMatriz A:\n{A}")
print(f"Tamaño de A: {A.shape}\n")

print(f"Vector B:\n{B}")
print(f"Tamaño de B: {B.shape}\n")

print(f"Solución U:\n{U}")
print(f"Tamaño de U: {U.shape}")

visualizar_resultado(nodes, elements, U)
