import numpy as np
from numpy.random import standard_normal
import matplotlib.pyplot as plt
from control.matlab import *  # MATLAB-like control toolbox functionality

# DC motor parameters
K = 4.9 # rad/s/V
tau = 0.085 # 1/s

# Process and measurement noise variance
q = 1000
r = 10
Q = q**2 * np.eye(2) #
#Q = 0 * np.eye(2)
R = np.eye(2) * r**2
#R = 0 * np.eye(2)

# Time between measurements, Sampling time, # iters.
Tm = 0.2 # Time between measurements
nm = 150  # number of measurements
ns = 5
Ts = Tm/ns # sampling time

# Initial conditions
x0  = [0,0] # true initial state
yhat0 = x0 + standard_normal() * np.sqrt(Q) # initial measurement
xhat0 = yhat0 # initial estimated state

# Setup KF system
A = np.array([[0, 1],
    [0, -1/tau]])
# A = np.array([0, 1],[0,-1/tau])
# print(A)
# A = np.array([0])
#A = -1/tau
B = np.array([[0], [K/tau]])
# B = K/tau
C = np.array([1, 0])
D = 0
sys_kf = ss(A, B, C, D)
# Setup covariance matrices
#Ad =  np.exp(A*Ts) #1 + A*Ts + A**2 * Ts**2/2
Ad =  np.eye(2) + A*Ts + A**2 * Ts**2/2
P0 = R

# Setup "true" system
Bn = np.hstack((B, np.eye(2)))
#print('bn', Bn)
Dn = np.array([0, 0, 0])
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
    zeta = np.matmul(standard_normal(size=(ns, 2)), np.sqrt(Q))
    zeta1 = zeta[:, 0]
    zeta2 = zeta[:, 1]
    un = np.stack((u, zeta1, zeta2),axis=1)

    y, _, x = lsim(sys_true, un, t, x0)
    #  add noise to sensor measurement
    yn = y[-1] + standard_normal() * np.sqrt(R)
    #  update initial condition
    x0 = x[-1]
    tn = T2
    # store values
    t_history = np.hstack((t_history, np.squeeze(t)))

    x_true_history = np.vstack((x_true_history, np.squeeze(x)))

    # Kalman filter
    #  Prediction
    y, _, x = lsim(sys_kf, u, t, xhat0[0])
    for i in range(ns-1):
        P = np.matmul(np.matmul(Ad, P), np.transpose(Ad)) + Ts**2 * Q
    # Store vals

    x_kf_history = np.vstack((np.squeeze(x_kf_history), np.squeeze(x)))

    #  Correction

    L = np.matmul(P, C) / (R + np.matmul(np.matmul(C, P), np.transpose(C)))
    x_plus = np.add(x[-1], np.matmul(L, (np.subtract(yn, y[-1]))))
    #print(np.shape)
    #P = (1 - L*C)*P*(1-L*C) + L*R*L
    P = np.add(np.matmul((np.matmul(np.subtract(np.eye(2), np.matmul(L, C)), P)), np.transpose(np.subtract(np.eye(2), np.matmul(L, C)))), np.matmul(np.matmul(L, R), np.transpose(L)))
    # Update IC
    xhat0 = x_plus

    print(np.shape(t_kf_plus))
    print(np.shape(T2))

"""    t_kf_plus = np.vstack((t_kf_plus, T2))
    x_kf_plus = np.vstack((x_kf_plus, x_plus))"""

fig, ax = plt.subplots()
x_kf_history = x_kf_history[:-1, :]

ax.plot(t_history, x_true_history[:, 0], label='true position')
ax.plot(t_history, x_kf_history[:, 0], label='KF position')
print(np.shape(t_kf_plus))
print(np.shape(x_kf_plus))
#ax.plot(t_kf_plus, x_kf_plus, 'rs',label='correction')
ax.set_xlabel('t [s]')
ax.set_ylabel(r'$\theta$ [rad]')
ax.legend()
plt.show()

fig, ax1 = plt.subplots()

ax1.plot(t_history, x_true_history[:, 1], label='true velocity')
ax1.plot(t_history, x_kf_history[:, 1], label='KF velocity')
#ax.plot(t_kf_plus, x_kf_plus, 'rs',label='correction')
ax1.set_xlabel('t [s]')
ax1.set_ylabel(r'$\omega$ [rad/s]')
ax1.legend()
plt.show()
