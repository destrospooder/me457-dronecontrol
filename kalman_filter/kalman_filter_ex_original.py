import numpy as np
from numpy.random import standard_normal
import matplotlib.pyplot as plt
from control.matlab import *  # MATLAB-like control toolbox functionality

# DC motor parameters 
K = 4.9 # rad/s/V 
tau = 0.085 # 1/s

# Process and measurement noise variance
Q  = 10**2 # 
R  = 0.01**2 # 1 standard dev.  

# Time between measurements, Sampling time, # iters.
Tm = 0.2 # Time between measurements
nm = 150  # number of measurements
ns = 5  
Ts = Tm/ns # sampling time

# Initial conditions
x0  = 0 # true initial state
yhat0  = x0 + standard_normal() * np.sqrt(Q) # initial measurement
xhat0 = yhat0 # initial estimated state

# Setup KF system
A = -1/tau
B = K/tau
C = 1
D = 0
sys_kf = ss(A, B, C, D)
# Setup covariance matrices
Ad =  np.exp(A*Ts) #1 + A*Ts + A**2 * Ts**2/2
P0 = R

# Setup "true" system
Bn = np.array([B, 1])
Dn = np.array([0, 0])
sys_true = ss(A, Bn, C, Dn)

# Init arrays for simulation
t_history = np.array([0])
x_true_history = np.array([x0])
x_kf_history = np.array([xhat0])
t_kf_plus = np.array([])
x_kf_plus = np.array([])
P_history = np.array([])
T1 = 0
T2 = 0
P = P0
for i in range(nm): # for all measurements
    T1 = T2
    T2 += Tm
    # Simulate "true" system for T1 to T2
    t = np.linspace(T1, T2, ns)
    u = np.sin(t) # np.ones(ns) 
    zeta = standard_normal(ns) * np.sqrt(Q)
    print(zeta)
    un = np.stack((u, zeta),axis=1)
    y, _, x = lsim(sys_true, un, t, x0)
    #  add noise to sensor measurement
    yn = y[-1] + standard_normal() * np.sqrt(R)
    #  update initial condition
    x0 = x[-1]
    tn = T2
    # store values
    t_history = np.hstack((t_history, np.squeeze(t)))
    x_true_history = np.hstack((x_true_history, np.squeeze(x)))
    # Kalman filter
    #  Prediction
    y, _, x = lsim(sys_kf, u, t, xhat0)
    for i in range(ns-1):
        P = Ad * P * Ad + Ts**2 * Q
    # Store vals
    x_kf_history = np.hstack((x_kf_history, np.squeeze(x)))
    #  Correction
    L = P * C / (R + C*P*C)
    x_plus = x[-1] + L*(yn - y[-1])
    P = (1 - L*C)*P*(1-L*C) + L*R*L
    # Update IC
    xhat0 = x_plus
    t_kf_plus = np.hstack((t_kf_plus, T2))
    x_kf_plus = np.hstack((x_kf_plus, x_plus))

fig, ax = plt.subplots()
ax.plot(t_history, x_true_history, label='true')
ax.plot(t_history, x_kf_history, label='KF')
ax.plot(t_kf_plus, x_kf_plus, 'rs',label='correction')
ax.set_xlabel('t [s]')
ax.set_ylabel(r'$\omega$ [rad/s]')
ax.legend()
plt.show()

