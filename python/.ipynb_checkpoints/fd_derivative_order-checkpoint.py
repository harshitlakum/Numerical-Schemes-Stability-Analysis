import numpy as np
import matplotlib.pyplot as plt

def u(x):
    return np.sin(x)

def u_prime_exact(x):
    return np.cos(x)

def fd_4th_order(x, h):
    return (-u(x + 2*h) + 8*u(x + h) - 8*u(x - h) + u(x - 2*h)) / (12*h)

x0 = 1.0
hs = 2.0 ** -np.arange(3, 11)
errors = np.abs(fd_4th_order(x0, hs) - u_prime_exact(x0))
orders = np.log(errors[:-1] / errors[1:]) / np.log(hs[:-1] / hs[1:])

print("   h         Error        Observed order")
for i in range(len(hs)):
    ordstr = f"{orders[i]:.3f}" if i < len(orders) else ""
    print(f"{hs[i]:.2e}   {errors[i]:.2e}   {ordstr}")

plt.figure()
plt.loglog(hs, errors, 'o-')
plt.xlabel('Step size h')
plt.ylabel('Absolute error')
plt.title('4th-Order FD Derivative Error vs h')
plt.grid(True, which='both')
plt.show()