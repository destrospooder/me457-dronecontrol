import numpy as np
import parameters as P
from integrators import get_integrator
from pid import PIDControl
import matplotlib.pyplot as plt


# f(t, x, u) = (Ku - y) / tau

def f(t, y, u):
    # mass sprint system state space form:
    return np.array([y[1], (1 / P.tau) * (P.K * u - y[1])])


class Controller:
    def __init__(self):
        kp = 5
        ki = 2.5
        kd = 1.5
        # TODO Look at saturation slope
        limit = P.umax
        sigma = P.sigma
        Ts = P.Ts
        self.controller = PIDControl(kp, ki, kd, limit, sigma, Ts, flag=False)

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

r = 1
y = np.array([0, 0])
t = 0

t_history = [0]
y_history = [y]
u_history = [0]

for i in range(P.nsteps):

    u = controller.update(r, y.item(0))
    y = system.update(u)
    t += P.Ts

    t_history.append(t)
    y_history.append(y)
    u_history.append(u)

# Plot response y due to step change in r
plt.close('all')

y_ = np.asarray(y_history)

fig, ax = plt.subplots()
ax.plot(t_history, y_[:,0], label='Motor Position', color='orange')
#ax.plot(t_history, y_[:,1], label='Motor Velocity', color='pink')
ax.set_xlabel('Time (s)')
ax.set_ylabel('Angular Position (rad)')
plt.axhline(y=r, color='r', linestyle='--', label='Commanded Position')
plt.legend()


ax2 = ax.twinx()
ax2.plot(t_history, u_history, label='Input Voltage', color='green')
ax2.set_ylabel('Actuation Signal (V)')
plt.legend(loc='lower right')
plt.show()

