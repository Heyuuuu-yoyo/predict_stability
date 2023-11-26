

import numpy as np
from scipy.io import loadmat
import matplotlib.pyplot as plt

# Load data
orange = loadmat('orange.mat')
purple = loadmat('purple.mat')
matlab1109 = loadmat('matlab1109.mat')

rr = matlab1109['r'].flatten()
A = matlab1109['A']
Day = 10
T = 200
# Initialization
S = 24
Time = 200
step = 0.1
Std_total = []
Stability = -np.ones((Time, 1))
Input_Abundance = np.zeros((Time, 54))
Input_RelativeAbundance = np.zeros((Time, 54))
Input_Presence = np.zeros((Time, 54))

# Simulation
for time in range(Time):
    Rand = np.random.permutation(54)
    species = np.sort(Rand[:S])
    AA = A[species, :][:, species]
    r = rr[species]

    N = np.zeros((Day*T+1, S))
    N[0, :] = np.random.rand(1, S)
    Input_Abundance[time, species] = N[0, :]
    Input_RelativeAbundance[time, species] = N[0, :] / np.sum(N[0, :])
    Input_Presence[time, species] = np.ones((1, S))

    for day in range(Day):
        for i in range(T*day + 1, T*day + T):
            for j in range(S):
                k1 = N[i-1, j] * (r[j] - np.dot(AA[j, :], N[i-1, :])) * step
                k2 = (N[i-1, j] + k1/2) * (r[j] - np.dot(AA[j, :], N[i-1, :])) * step
                k3 = (N[i-1, j] + k2/2) * (r[j] - np.dot(AA[j, :], N[i-1, :])) * step
                k4 = (N[i-1, j] + k3) * (r[j] - np.dot(AA[j, :], N[i-1, :])) * step
                N[i, j] = N[i-1, j] + (1/6) * (k1 + 2*k2 + 2*k3 + k4) + 10**-6 * step
                if N[i, j] > 1:
                    N[i, j] = 1

    Std = np.std(N[1800:2000, :], axis=0)
    Std_total.append(np.max(Std))
    if (np.max(Std) > 0.01):
        Stability[time] = 1
    else:
        Stability[time] = 0

print(Stability)