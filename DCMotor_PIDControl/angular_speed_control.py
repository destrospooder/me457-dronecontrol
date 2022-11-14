import numpy as np
import parameters as P
from integrators import get_integrator
from pid import PIDControl
import matplotlib.pyplot as plt


# f(t, x, u) = (Ku - y) / tau

def f(t, y, u):
    return (P.K * u - y) / P.tau


class Controller:
    def __init__(self):
        kp = P.kp
        ki = P.ki
        kd = P.kd
        # TODO Look at saturation slope
        limit = P.emax
        sigma = P.sigma
        Ts = P.Ts
        self.controller = PIDControl(kp, ki, kd, limit, sigma, Ts)

    def update(self, r, y):
        return self.controller.PID(r, y)


class System:
    def __init__(self):
        self.integrator = get_integrator(P.Ts, f)
        self.t_history = [0]
        self.y_history = [0]
        pass

    def update(self, u):
        return self.integrator.step(t, y, u)

        pass


# Init system and feedback controller
system = System()
controller = Controller()

# Simulate step response

r = 5
t = 0
y = 0
u = 0


print(P.kp)
print(P.ki)

t_history = [t]
y_history = [y]
u_history = [u]

for i in range(P.nsteps):
    u = controller.update(r, y)
    y = system.update(u)
    t += P.Ts

    t_history.append(t)
    y_history.append(y)
    u_history.append(u)

# Plot response y due to step change in r
plt.close('all')


fig, ax = plt.subplots()
ax.plot(t_history, y_history, label='Motor Speed', color='orange')
ax.set_xlabel('Time (s)')
ax.set_ylabel('Angular Velocity (rad/s)')
plt.axhline(y=r, color='r', linestyle='--', label='Commanded Speed')
plt.legend()


ax2 = ax.twinx()
ax2.plot(t_history, u_history, label='Input Voltage', color='green')
ax2.set_ylabel('Actuation Signal (V)')
plt.legend()
plt.show()

