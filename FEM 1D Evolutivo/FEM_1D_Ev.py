import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

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
    return np.cos(2*np.pi*x)

# ---------------------
# Malla
# ---------------------
x = np.linspace(0, L, Nnodes)

# Matrices globales
K = np.zeros((Nnodes, Nnodes))
M = np.zeros((Nnodes, Nnodes))
F = np.zeros(Nnodes)

# ---------------------
# Cuadratura Gauss
# ---------------------
xi = np.array([-1/np.sqrt(3), 1/np.sqrt(3)])
w = np.array([1.0, 1.0])

# ---------------------
# Ensamblaje FEM
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
# Dirichlet
# ---------------------
def apply_dirichlet(A, b):
    # x=0
    A[0,:] = 0; A[:,0] = 0
    A[0,0] = 1; b[0] = 0

    # x=L
    A[-1,:] = 0; A[:,-1] = 0
    A[-1,-1] = 1; b[-1] = 1

# ---------------------
# Condición inicial
# ---------------------
U = np.zeros(Nnodes)
U[-1] = 1.0 # Aplicamos Dirichlet desde el inicio

U_history = [U.copy()] # Guardamos el estado inicial t=0

# Pre-calculamos A si no cambia en el tiempo para ganar velocidad
A_base = M + dt * K 

for n in range(Nt):
    # b depende de la U del paso anterior
    b = M @ U + dt * F 
    
    # Creamos una copia de A para no destruir la original con Dirichlet
    A_step = A_base.copy()
    
    apply_dirichlet(A_step, b)

    U = np.linalg.solve(A_step, b)
    U_history.append(U.copy())

# Nota: Ahora U_array tendrá tamaño (Nt + 1, Nnodes)
U_array = np.array(U_history)

# ---------------------
# Bucle temporal
# ---------------------
U_history = []

for n in range(Nt):

    A = M + dt*K
    b = M @ U + dt*F

    apply_dirichlet(A, b)

    U = np.linalg.solve(A, b)

    U_history.append(U.copy())

U_array = np.array(U_history)

print(f"Valor del vector solución en t={Tfinal:.2f}: {U}")
print(f"Tamaño de U_array: {U_array.shape}")

# =====================
# GRÁFICAS
# =====================

# ---------------------
# Curvas en distintos tiempos
# ---------------------
plt.figure()

times_to_plot = [0, Nt//4, Nt//2, Nt-1]

for t in times_to_plot:
    plt.plot(x, U_array[t], label=f"t={t*dt:.2f}")

plt.xlabel("x")
plt.ylabel("u(x,t)")
plt.title("Evolución temporal (curvas)")
plt.legend()
plt.grid()

# Guardar figura
plt.savefig("curvas.png")
plt.close()

# ---------------------
# Mapa de calor
# ---------------------
plt.figure()

plt.imshow(U_array,
           aspect='auto',
           extent=[x[0], x[-1], 0, Tfinal],
           origin='lower')

plt.colorbar(label="u(x,t)")
plt.xlabel("x")
plt.ylabel("t")
plt.title("Mapa de calor")

# Guardar figura
plt.savefig("heatmap.png")
plt.close()

# ---------------------
# Animación
# ---------------------
fig, ax = plt.subplots()
line, = ax.plot(x, U_array[0])

ax.set_xlim(x[0], x[-1])
ax.set_ylim(np.min(U_array), np.max(U_array))

def update(frame):
    line.set_ydata(U_array[frame])
    ax.set_title(f"t = {frame*dt:.2f}")
    return line,

ani = FuncAnimation(fig, update, frames=Nt, interval=50)

# Guardar animación (GIF)
ani.save("animacion.gif", writer="pillow", fps=30)

print("Animación guardada como animacion.gif")
