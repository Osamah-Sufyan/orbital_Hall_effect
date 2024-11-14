

import numpy as np
import matplotlib.pyplot as plt

######## orbital Berry curvature ########

m_e = 9.109 *10**-31 # electron mass
h_bar = 1.054571 *10**-34
g_L = 2.002319 #Lande g-factor
a = 0.246* 10**-9 #lattice constant
delta = 1*10**-19
t = 4.325*10**-19

def f_k(kx, ky):
    term1 = np.exp(-1j * ky * a)
    term2 = 1 + 2 * np.cos(np.sqrt(3) * kx * a / 2) * np.exp(1j * 3 * ky * a / 2)
    return t * term1 * term2

def omega_orbit(kx, ky, a, t):


    prefactor = -2 * m_e / (8 * g_L * np.square(h_bar))
    
    term1 = np.square(3 * np.sqrt(3) * np.square(a)*np.square(t) * delta) / 16
    
    term2 = np.square(np.sin(np.sqrt(3) * kx * a))
    
    denominator = np.power(np.square(delta) + np.square(np.abs(f_k(kx,ky))), 5/2)
    
    return prefactor * term1 * term2 / denominator

# kx = np.linspace(-np.pi / a, np.pi / a, 100)
# ky = np.linspace(-np.pi / a, np.pi / a, 100)
# kx, ky = np.meshgrid(kx, ky)
# ok = omega_orbit(kx, ky, a, t)
# plt.imshow(np.abs(ok), extent=(-np.pi / a, np.pi / a, -np.pi / a, np.pi / a))
# plt.colorbar()
# plt.xlabel('kx')
# plt.ylabel('ky')
# plt.title('f(k)')
# plt.show()



######## Berry curvature ########
def omega(kx, ky, a, t):


    prefactor = -2 * m_e / (8 * g_L * np.square(h_bar))
    
    term1 = 3 * np.sqrt(3) * np.square(a)*np.square(t) * delta / 4
    
    term2 = np.sin(np.sqrt(3) * kx * a)
    
    denominator = np.power(np.square(delta) + np.square(np.abs(f_k(kx,ky))), 3/2)
    
    return  term1 * term2 / denominator

kx = np.linspace(-np.pi / a, np.pi / a, 100)
ky = np.linspace(-np.pi / a, np.pi / a, 100)
kx, ky = np.meshgrid(kx, ky)
ok = omega(kx, ky, a, t)
plt.imshow(ok, extent=(-np.pi / a, np.pi / a, -np.pi / a, np.pi / a))
plt.colorbar()
plt.xlabel(r'$k_x$', fontsize=16)
plt.ylabel(r'$k_y$', fontsize=16)
plt.title('normal Berry Curvature')
plt.savefig('norml_berry_curvature.png')
plt.show()