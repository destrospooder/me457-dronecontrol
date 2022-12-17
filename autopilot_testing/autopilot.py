import sys
import control
import numpy as np
from numpy import array, sin, cos, radians, concatenate, zeros, diag
from scipy.linalg import solve_continuous_are, inv

import control_parameters as AP
from wrap import wrap
import model_coef as M
from msg_state import MsgState
from msg_delta import MsgDelta
import matplotlib.pyplot as plt
sys.path.append('..')

class Autopilot:
    def __init__(self, ts_control):
        self.Ts = ts_control

        self.integratorCourse = 0
        self.integratorAltitude = 0
        self.integratorAirspeed = 0
        self.errorCourseD1 = 0
        self.errorAltitudeD1 = 0
        self.errorAirspeedD1 = 0

        CrLat = array([[0, 0, 0, 0, 1.0]])
        AAlat = concatenate((
                    concatenate((M.A_lat, zeros((5, 1))), axis=1),
                    concatenate((CrLat, zeros((1, 1))), axis=1)),
                axis=0)
        BBlat = concatenate((M.B_lat, zeros((1, 2))), axis=0)
        Qlat = diag([0.001, 0.01, 0.1, 100, 1, 100]) # v, p, r, phi, chi, intChi
        Rlat = diag([1, 1]) #a, r
        Plat = solve_continuous_are(AAlat, BBlat, Qlat, Rlat)
        self.Klat = inv(Rlat) @ BBlat.T @ Plat
        CrLon = array([[0, 0, 0, 0, 1.0], [1/AP.Va0, 1/AP.Va0, 0, 0, 0]])
        AAlon = concatenate((
                    concatenate((M.A_lon, zeros((5, 2))), axis=1),
                    concatenate((CrLon, zeros((2, 2))), axis=1)),
                axis=0)
        BBlon = concatenate((M.B_lon, zeros((2, 2))), axis=0)
        Qlon = diag([10, 10, 0.001, 0.01, 10, 100, 100])
        Rlon = diag([1, 1])
        Plon = solve_continuous_are(AAlon, BBlon, Qlon, Rlon)
        self.Klon = inv(Rlon) @ BBlon.T @ Plon
        self.commanded_state = MsgState

    def update(self, cmd, state):
        # lateral autopilot
        errorAirspeed = state.Va - cmd.airspeed_command
        chi_c = wrap(cmd.course_command, state.chi)
        errorCourse = saturate(state.chi - chi_c, -radians(15), radians(15))
        self.integratorCourse = self.integratorCourse + (self.Ts/2) * (errorCourse + self.errorCourseD1)
        self.errorCourseD1 = errorCourse
        xLat = array([[errorAirspeed * sin(state.beta)],
                      [state.p],
                      [state.r],
                      [state.phi],
                      [errorCourse],
                      [self.integratorCourse]])
        tmp = -self.Klat @ xLat
        delta_a = saturate(tmp.item(0), -radians(30), radians(30))
        delta_r = saturate(tmp.item(1), -radians(30), radians(30))

        # long. autopilot
        altitude_c = saturate(cmd.altitude_command,
                              state.altitude - 0.2*AP.altitude_zone,
                              state.altitude + 0.2*AP.altitude_zone)

        errorAltitude = state.altitude - altitude_c
        self.integratorAltitude = self.integratorAltitude + (self.Ts/2) * (errorAltitude + self.errorAltitudeD1)
        self.errorAltitudeD1 = errorAltitude
        self.integratorAirspeed = self.integratorAirspeed + (self.Ts / 2) * (errorAirspeed + self.errorAirspeedD1)
        self.integratorAirspeed = errorAirspeed
        xLon = array([[errorAirspeed * cos(state.alpha)],
                      [errorAirspeed * sin(state.alpha)],
                      [state.q],
                      [state.theta],
                      [errorAltitude],
                      [self.integratorAltitude],
                      [self.integratorAirspeed]])

        tmp = -self.Klon @ xLon
        delta_e = saturate(tmp.item(0), -radians(30), radians(30))
        delta_t = saturate(tmp.item(1), 0.0, 1.0)

        delta = MsgDelta(elevator=delta_e,
                         aileron=delta_a,
                         rudder=delta_r,
                         throttle=delta_t)
        self.commanded_state.altitude = cmd.altitude_command
        self.commanded_state.Va = cmd.airspeed_command
        return delta, self.commanded_state

    def step_lat(self):

        CrLat = array([[0, 0, 0, 0, 1.0]])
        AAlat = concatenate((
            concatenate((M.A_lat, zeros((5, 1))), axis=1),
            concatenate((CrLat, zeros((1, 1))), axis=1)),
            axis=0)
        BBlat = concatenate((M.B_lat, zeros((1, 2))), axis=0)
        Qlat = diag([0.001, 0.01, 0.1, 100, 1, 100])  # v, p, r, phi, chi, intChi
        Rlat = diag([1, 1])  # a, r
        Plat = solve_continuous_are(AAlat, BBlat, Qlat, Rlat)
        self.Klat = inv(Rlat) @ BBlat.T @ Plat

        system = control.ss(np.subtract(AAlat, np.matmul(BBlat, self.Klat)), BBlat, CrLat, 0)

        y, t = control.step_response(system)

        plt.close('all')

        fig, ax = plt.subplots()
        ax.plot(t, y, label='Motor Position', color='orange')
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Angular Position (rad)')

        plt.show()





def saturate(input, low_limit, up_limit):
    if input <= low_limit:
        output = low_limit
    elif input >= up_limit:
        output = up_limit
    else:
        output = input
    return output


"""def df_dx(x, z):
    eps = 0.001
    A = np.zeros((m, n))
    f_at_x = f(x, z)
    for i in range(0, n):
        x_eps = np.copy(x)
        x_eps[i][0] += eps
        f_at_x_eps = f(x_eps, z)
        df_dxi = (f_at_x_eps - f_at_x) / eps
        A[;, i] = df_dxi[:, 0]

    return A
"""