import numpy as np

def stable_matrix(m, lam=0.9, seed=0):
    rng = np.random.default_rng(seed)
    Q, _ = np.linalg.qr(rng.standard_normal((m, m)))
    diag_vals = lam * rng.uniform(0.2, 1.0, size=m)
    D = Q @ np.diag(diag_vals) @ np.linalg.inv(Q)
    return D

def run_scheme(D, f_list, phi):
    u = phi.copy()
    us = [u.copy()]
    for f in f_list:
        u = D @ u + f
        us.append(u.copy())
    return us

def power_norm_bound(D, nmax):
    m = D.shape[0]
    C = 0.0
    P = np.eye(m)
    for k in range(nmax+1):
        C = max(C, np.linalg.norm(P, 2))
        P = D @ P
    return C

def main():
    m = 40
    nsteps = 100
    D = stable_matrix(m, lam=0.95, seed=42)
    rng = np.random.default_rng(123)
    phi = rng.standard_normal(m)
    f_list = [0.01 * rng.standard_normal(m) for _ in range(nsteps)]
    us = run_scheme(D, f_list, phi)
    CT = power_norm_bound(D, nsteps)
    sum_f = np.sum([np.linalg.norm(f, 2) for f in f_list])
    rhs = CT * (np.linalg.norm(phi, 2) + sum_f)
    left = np.array([np.linalg.norm(u, 2) for u in us])
    holds = np.all(left <= rhs + 1e-12)
    print(f"C_T (max power norm): {CT:.3e}")
    print(f"||phi||_2 = {np.linalg.norm(phi):.3e}")
    print(f"sum ||f_j||_2 = {sum_f:.3e}")
    print(f"RHS bound = {rhs:.3e}")
    print(f"Max ||u^n||_2 = {left.max():.3e}")
    print("Inequality holds for all n:", holds)
    for n in [0, nsteps//2, nsteps]:
        print(f"n={n:3d}: ||u^n|| = {left[n]:.3e}  <=  RHS = {rhs:.3e}")

if __name__ == "__main__":
    main()