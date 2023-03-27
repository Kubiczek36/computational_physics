# expected height of a particle in homogeneous gravitational field calculated by Monte Carlo method

# %% import modules
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import Boltzmann as kB
from scipy.constants import g

# %% define constants
m = 4.66e-26 # mass of the particle [kg]
T = 300.0 # temperature [K]
beta = 1/(kB*T) # inverse temperature [1/K]

def f(h):
    return np.exp(-beta*g*h)

# %% Monte Carlo method

N = 1 # number of particles
h = np.linspace(0, stop=1e-8, num=1000) # height of the particle [m]

# %% calculate the expected height of the particle

def rho(h):
    return np.exp(-beta*g*h)

# plot the probability density function
plt.plot(h, rho(h))
plt.xlabel('h [m]')
plt.ylabel(r'$\rho (h)$')
plt.show()

# %%

exp_g = 1/N * np.sum(rho(h))
print(f"expected height is {exp_g}")

# try sampling with np.random.exp
# https://numpy.org/doc/stable/reference/random/generated/numpy.random.exponential.html