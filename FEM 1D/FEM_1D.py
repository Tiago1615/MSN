import numpy as np

# Parámetros del problema
L = 1.0
Nnodes = 5
c = 1.0
u0 = 0.0
uL_prime = 0.0

def f(x):
    return 1.0

# ---------------------
# Malla
# ---------------------
x = np.linspace(0, L, Nnodes)
Nelements = Nnodes - 1

# Conectividad
elem = np.zeros((Nelements,2), dtype=int)
for n in range(Nelements):
    elem[n,0] = n
    elem[n,1] = n+1

# ---------------------
# Inicializar sistema global
# ---------------------
A = np.zeros((Nnodes, Nnodes))
B = np.zeros(Nnodes)

# ---------------------
# Cuadratura gaussiana 2 puntos
# ---------------------
xi = np.array([(1-1/np.sqrt(3))/2,
               (1+1/np.sqrt(3))/2])
w = np.array([0.5, 0.5])

# ---------------------
# Bucle en elementos
# ---------------------
for n in range(Nelements):
    
    nodes = elem[n]
    x1 = x[nodes[0]]
    x2 = x[nodes[1]]
    h = x2 - x1
    
    # Inicializar matrices elementales
    Ae = np.zeros((2,2))
    Be = np.zeros(2)
    
    # Bucle en puntos de integración
    for k in range(2):
        
        xi_k = xi[k]
        w_k = w[k]
        
        # Funciones de forma en referencia
        Nhat = np.array([1-xi_k, xi_k])
        
        # Derivadas en referencia
        dNhat = np.array([-1.0, 1.0])
        
        # Transformación al elemento físico
        x_phys = x1 + h*xi_k
        
        # Derivadas físicas
        dNdx = dNhat / h
        
        # Construcción matriz elemental
        for alpha in range(2):
            for beta in range(2):
                
                # término rigidez
                Ae[alpha,beta] += w_k * (
                    (1/h)*dNhat[alpha]*dNhat[beta]
                )
                
                # término masa
                Ae[alpha,beta] += w_k * (
                    c*h*Nhat[alpha]*Nhat[beta]
                )
            
            # Vector elemental
            Be[alpha] += w_k * (
                h * f(x_phys) * Nhat[alpha]
            )
    
    # -----------------
    # Ensamblaje
    # -----------------
    for alpha in range(2):
        I = nodes[alpha]
        B[I] += Be[alpha]
        
        for beta in range(2):
            J = nodes[beta]
            A[I,J] += Ae[alpha,beta]

# ---------------------
# Condición Neumann
# ---------------------
B[-1] += uL_prime

# ---------------------
# Dirichlet en nodo 0
# ---------------------
A[0,:] = 0
A[0,0] = 1
B[0] = u0

# ---------------------
# Resolver sistema
# ---------------------
U = np.linalg.solve(A,B)

print(U)