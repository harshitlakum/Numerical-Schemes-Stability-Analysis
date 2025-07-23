import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags, kron, identity
from scipy.sparse.linalg import spsolve

def u_exact(x, y):
    return np.sin(x * y)

def f_rhs(x, y):
    return (x**2 + y**2) * np.sin(x * y)

def solve_poisson_dirichlet(f, g, Nx, Ny=None):
    if Ny is None:
        Ny = Nx
    hx = 1.0 / (Nx + 1)
    hy = 1.0 / (Ny + 1)
    x = np.linspace(hx, 1.0 - hx, Nx)
    y = np.linspace(hy, 1.0 - hy, Ny)
    Xint, Yint = np.meshgrid(x, y, indexing='ij')
    Tx = diags([-1, 2, -1], offsets=[-1,0,1], shape=(Nx, Nx)) / (hx*hx)
    Ty = diags([-1, 2, -1], offsets=[-1,0,1], shape=(Ny, Ny)) / (hy*hy)
    Ix = identity(Nx)
    Iy = identity(Ny)
    A = -(kron(Iy, Tx) + kron(Ty, Ix))
    F = f(Xint, Yint).ravel()
    xb = np.linspace(0, 1, Nx+2)
    yb = np.linspace(0, 1, Ny+2)
    for i in range(Nx):
        xi = xb[i+1]
        F[i*Ny + 0]   -= g(xi, 0.0)   / (hy*hy)
        F[i*Ny + (Ny-1)] -= g(xi, 1.0) / (hy*hy)
    for j in range(Ny):
        yj = yb[j+1]
        F[0*Ny + j]   -= g(0.0, yj)   / (hx*hx)
        F[(Nx-1)*Ny + j] -= g(1.0, yj) / (hx*hx)
    U_int = spsolve(A, F)
    U = np.zeros((Nx+2, Ny+2))
    U[1:-1, 1:-1] = U_int.reshape((Nx, Ny), order='C')
    Xfull, Yfull = np.meshgrid(xb, yb, indexing='ij')
    U[0,  :] = g(0.0,  yb)
    U[-1, :] = g(1.0,  yb)
    U[:,  0] = g(xb,  0.0)
    U[:, -1] = g(xb,  1.0)
    return Xfull, Yfull, U

def plot_solution_and_error(X, Y, Uh):
    Uex = u_exact(X, Y)
    Err = Uh - Uex
    fig1 = plt.figure(figsize=(6,5))
    ax1 = fig1.add_subplot(111, projection='3d')
    ax1.plot_surface(X, Y, Uh, linewidth=0, antialiased=True)
    ax1.set_title("Numerical solution $u_h$")
    ax1.set_xlabel("x"); ax1.set_ylabel("y")
    fig2 = plt.figure(figsize=(6,5))
    ax2 = fig2.add_subplot(111, projection='3d')
    ax2.plot_surface(X, Y, np.abs(Err), linewidth=0, antialiased=True)
    ax2.set_title("Pointwise error $|u_h - u|$")
    ax2.set_xlabel("x"); ax2.set_ylabel("y")
    plt.show()

def convergence_test():
    Ns = [16, 32, 64, 128]
    errs_inf = []
    errs_L2  = []
    for N in Ns:
        X, Y, Uh = solve_poisson_dirichlet(f_rhs, u_exact, N)
        Uex = u_exact(X, Y)
        err = np.abs(Uh - Uex)
        errs_inf.append(err.max())
        hx = 1/(N+1); hy = 1/(N+1)
        errs_L2.append(np.sqrt(np.sum(err[1:-1,1:-1]**2)*hx*hy))
    def orders(errs):
        errs = np.array(errs, dtype=float)
        return np.log(errs[:-1]/errs[1:]) / np.log(np.array(Ns[:-1])/np.array(Ns[1:]))
    ord_inf = orders(errs_inf)
    ord_L2  = orders(errs_L2)
    print("N    ||e||_inf        ||e||_2")
    for N, ei, e2 in zip(Ns, errs_inf, errs_L2):
        print(f"{N:4d}  {ei:12.4e}  {e2:12.4e}")
    print("\nObserved order (inf-norm):", ord_inf)
    print("Observed order (L2-norm) :", ord_L2)
    return Ns, errs_inf, errs_L2

def plot_convergence(Ns, errs_inf, errs_L2):
    Ns = np.array(Ns, dtype=float)
    plt.figure()
    plt.loglog(1.0/Ns, errs_inf, 'o-', label='inf-norm')
    plt.loglog(1.0/Ns, errs_L2,  's-', label='L2-norm')
    plt.xlabel('h = 1/(N+1)')
    plt.ylabel('error')
    plt.title('Convergence of 2D FD Poisson solver')
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.legend()
    plt.show()

if __name__ == "__main__":
    X, Y, Uh = solve_poisson_dirichlet(f_rhs, u_exact, Nx=64)
    plot_solution_and_error(X, Y, Uh)
    Ns, e_inf, e_L2 = convergence_test()
    plot_convergence(Ns, e_inf, e_L2)