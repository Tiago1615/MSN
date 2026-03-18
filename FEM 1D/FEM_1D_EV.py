import numpy as np

# ---------------------
# Parámetros
# ---------------------
L = 10.0
Nelements = 100
Nnodes = Nelements + 1

alpha = 1.0
dt = 0.01
Tfinal = 1.0
Nt = int(Tfinal/dt)

def f(x):
    return 2.0 * x

# ---------------------
# Malla
# ---------------------
x = np.linspace(0, L, Nnodes)

K = np.zeros((Nnodes, Nnodes))
M = np.zeros((Nnodes, Nnodes))
F = np.zeros(Nnodes)

# ---------------------
# Cuadratura Gauss
# ---------------------
xi = np.array([-1/np.sqrt(3), 1/np.sqrt(3)])
w = np.array([1.0, 1.0])

# ---------------------
# Ensamblaje
# ---------------------
for n in range(Nelements):

    x1 = x[n]
    x2 = x[n+1]
    h = x2 - x1

    Ke = np.zeros((2,2))
    Me = np.zeros((2,2))
    Fe = np.zeros(2)

    for k in range(2):

        xi_k = xi[k]
        w_k = w[k]

        Nhat = np.array([(1 - xi_k)/2,
                         (1 + xi_k)/2])
        dNhat = np.array([-0.5, 0.5])

        x_phys = (x1 + x2)/2 + (h/2)*xi_k
        J = h/2
        dNdx = dNhat * (2/h)

        for i in range(2):
            for j in range(2):
                Ke[i,j] += w_k * dNdx[i]*dNdx[j]*J
                Me[i,j] += w_k * Nhat[i]*Nhat[j]*J

            Fe[i] += w_k * f(x_phys)*Nhat[i]*J

    # ensamblaje global
    nodes = [n, n+1]
    for i in range(2):
        I = nodes[i]
        F[I] += Fe[i]

        for j in range(2):
            Jg = nodes[j]
            K[I,Jg] += Ke[i,j]
            M[I,Jg] += Me[i,j]

# ---------------------
# Condiciones iniciales
# ---------------------
U = np.zeros(Nnodes)  # u(x,0)

# ---------------------
# Condiciones Dirichlet
# ---------------------
def apply_dirichlet(A, b):
    A[0,:] = 0; A[:,0] = 0
    A[0,0] = 1; b[0] = 0

    A[-1,:] = 0; A[:,-1] = 0
    A[-1,-1] = 1; b[-1] = 0

# ---------------------
# Bucle temporal
# ---------------------
for n in range(Nt):

    A = M + dt*K
    b = M @ U + dt*F

    apply_dirichlet(A, b)

    U = np.linalg.solve(A, b)

print(U)