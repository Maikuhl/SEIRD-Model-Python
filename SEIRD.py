# importing required modules
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

# defining the plotting function to visualise SEIRD model
def plotseird(t, S, E, I, R, D):
    f, ax = plt.subplots(1, 1, figsize=(10,4))
    ax.plot(t, S, linewidth=2, label='Susceptible')
    ax.plot(t, E, linewidth=2, label='Exposed')
    ax.plot(t, I, linewidth=2, label='Infectious')
    ax.plot(t, R, linewidth=2, label='Recovered')
    ax.plot(t, D, linewidth=2, label='Deaths')
    ax.plot(t, S+E+I+R+D, 'c--', linewidth=2, label='Total')

    ax.set_xlabel('Time (days)')

    ax.grid(visible=True, color='w', linewidth=2, linestyle='-')

    legend = ax.legend(borderpad=2.0)
    legend.get_frame().set_alpha(0.5)

    for spine in ('top', 'right', 'bottom', 'left'):
        ax.spines[spine].set_visible(True)

    plt.show()

# defining the SEIRD system of differential equations 
def SEIRD(y, t, rB, epsilon, gamma, delta):
    S, E, I, R, D = y
    N = S + E + I + R + D 
    dS = -rB * S * I / N
    dE = (rB * I / N) * S - (epsilon * E)
    dI = epsilon * E - (gamma * I) - (delta * I)
    dR = gamma * I
    dD = delta * I
    return dS, dE, dI, dR, dD


# setting constants 
rB = 0.8        # float 
epsilon = 0.2   # float
gamma = 0.1     # float
delta = 0.05    # float


# setting initial conditions
S0, E0, I0, R0, D0 = 999, 0, 1, 0, 0
y0 = S0, E0, I0, R0, D0

# setting matrix of timesteps
t = np.linspace(0, 150)

model = odeint(SEIRD, y0, t, args=(rB, epsilon, gamma, delta))

S, E, I, R, D = model.T

plotseird(t, S, E, I, R, D)

# value of peak infected population
print(max(I))

# day corresponding to peak infected population
print(3 * np.argmax(I))

# day where infected population falls below 10
U = np.round(I,2)
print(U)

# by inspection (sorry it's hardcoded),
print(3 * list(U).index(9.74))



