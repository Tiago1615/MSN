import numpy as np

# ---------------------
# Parámetros
# ---------------------
L = 10.0
Nelements = 100
Nnodes = Nelements + 1  # 101 nodos

c = 1.0
u0 = 0.0
uL_prime = 1.0

def f(x):
    return 2.0 * x

# ---------------------
# Malla
# ---------------------
x = np.linspace(0, L, Nnodes)

# CAMBIO: Dimensiones basadas en Nnodes
A = np.zeros((Nnodes, Nnodes))
B = np.zeros(Nnodes)

# ---------------------
# Cuadratura Gauss [-1,1]
# ---------------------
xi = np.array([-1/np.sqrt(3), 1/np.sqrt(3)])
w = np.array([1.0, 1.0])

# ---------------------
# Bucle en elementos
# ---------------------
# CAMBIO: range(Nelements) para incluir el último elemento
for n in range(Nelements):

    x1 = x[n]
    x2 = x[n+1]
    h = x2 - x1

    Ae = np.zeros((2,2))
    Be = np.zeros(2)

    for k in range(2):
        xi_k = xi[k]
        w_k = w[k]

        Nhat = np.array([(1 - xi_k)/2, (1 + xi_k)/2])
        dNhat = np.array([-0.5, 0.5])

        x_phys = (x1 + x2)/2 + (h/2)*xi_k
        J = h/2
        dNdx = dNhat * (1/J) # Es lo mismo que dNhat * (2/h)

        for alpha in range(2):
            for beta in range(2):
                # Rigidez + Masa
                Ae[alpha,beta] += w_k * (dNdx[alpha]*dNdx[beta]*J)
                Ae[alpha,beta] += w_k * (c*Nhat[alpha]*Nhat[beta]*J)

            Be[alpha] += w_k * (f(x_phys)*Nhat[alpha]*J)

    # -----------------
    # Ensamblaje
    # -----------------
    global_nodes = [n, n+1]

    for alpha in range(2):
        I = global_nodes[alpha]
        B[I] += Be[alpha]
        for beta in range(2):
            Jg = global_nodes[beta]
            # CAMBIO: Ya no necesitas el "if Jg < Nelements"
            A[I,Jg] += Ae[alpha,beta]

# ---------------------
# Neumann (en el último nodo: Nnodes - 1)
# ---------------------
B[-1] += uL_prime

# ---------------------
# Dirichlet (en el primer nodo: 0)
# ---------------------
A[0, :] = 0.0
A[0, 0] = 1.0
B[0] = u0
# Nota: La eliminación de simetría (A[:, 0] = 0) es opcional si fijas bien la fila.

# ---------------------
# Resolver
# ---------------------
U = np.linalg.solve(A,B)
U = U + 0.0  # Para evitar -0.0

np.set_printoptions(precision=6, suppress=True)
print(np.round(U, 10))