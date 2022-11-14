"""
mavDynamics 
    - this file implements the dynamic equations of motion for MAV
    - use unit quaternion for the attitude state
    
part of mavPySim 
    - Beard & McLain, PUP, 2012
    - Update history:  
        12/20/2018 - RWB
"""
import sys
sys.path.append('..')
import numpy as np
import math
# load message types
from message_types.msg_state import MsgState

import parameters.aerosonde_parameters as MAV
from tools.rotations import Quaternion2Rotation, Quaternion2Euler


class MavDynamics:
    def __init__(self, Ts):
        self._ts_simulation = Ts
        # set initial states based on parameter file
        # _state is the 13x1 internal state of the aircraft that is being propagated:
        # _state = [pn, pe, pd, u, v, w, e0, e1, e2, e3, p, q, r]
        # We will also need a variety of other elements that are functions of the _state and the wind.
        # self.true_state is a 19x1 vector that is estimated and used by the autopilot to control the aircraft:
        # true_state = [pn, pe, h, Va, alpha, beta, phi, theta, chi, p, q, r, Vg, wn, we, psi, gyro_bx, gyro_by, gyro_bz]
        self._state = np.array([[MAV.north0],  # (0)
                               [MAV.east0],   # (1)
                               [MAV.down0],   # (2)
                               [MAV.u0],    # (3)
                               [MAV.v0],    # (4)
                               [MAV.w0],    # (5)
                               [MAV.e0],    # (6)
                               [MAV.e1],    # (7)
                               [MAV.e2],    # (8)
                               [MAV.e3],    # (9)
                               [MAV.p0],    # (10)
                               [MAV.q0],    # (11)
                               [MAV.r0]])   # (12)
        # store wind data for fast recall since it is used at various points in simulation
        self._wind = np.array([[0.], [0.], [0.]])  # wind in NED frame in meters/sec
        self._update_velocity_data()
        # store forces to avoid recalculation in the sensors function
        self._forces = np.array([[0.], [0.], [0.]])
        self._Va = MAV.u0
        self._alpha = 0
        self._beta = 0
        # initialize true_state message
        self.true_state = MsgState()

    ###################################
    # public functions
    def update(self, delta, wind):
        """
            Integrate the differential equations defining dynamics, update sensors
            delta = (delta_a, delta_e, delta_r, delta_t) are the control inputs
            wind is the wind vector in inertial coordinates
            Ts is the time step between function calls.
        """
        # get forces and moments acting on rigid bod
        forces_moments = self._forces_moments(delta)

        # Integrate ODE using Runge-Kutta RK4 algorithm
        time_step = self._ts_simulation
        k1 = self._derivatives(self._state, forces_moments)
        k2 = self._derivatives(self._state + time_step/2.*k1, forces_moments)
        k3 = self._derivatives(self._state + time_step/2.*k2, forces_moments)
        k4 = self._derivatives(self._state + time_step*k3, forces_moments)
        self._state += time_step/6 * (k1 + 2*k2 + 2*k3 + k4)

        # normalize the quaternion
        e0 = self._state.item(6)
        e1 = self._state.item(7)
        e2 = self._state.item(8)
        e3 = self._state.item(9)
        normE = (e0**2+e1**2+e2**2+e3**2) ** (1/2)
        self._state[6][0] = self._state.item(6)/normE
        self._state[7][0] = self._state.item(7)/normE
        self._state[8][0] = self._state.item(8)/normE
        self._state[9][0] = self._state.item(9)/normE

        # update the airspeed, angle of attack, and sideslip angles using new state
        self._update_velocity_data(wind)

        # update the message class for the true state
        self._update_true_state()

    def external_set_state(self, new_state):
        self._state = new_state

    ###################################
    # private functions
    def _derivatives(self, state, forces_moments):
        """
        for the dynamics xdot = f(x, u), returns f(x, u)
        """
        # extract the states
        north = state.item(0)
        east = state.item(1)
        down = state.item(2)
        u = state.item(3)
        v = state.item(4)
        w = state.item(5)
        e0 = state.item(6)
        e1 = state.item(7)
        e2 = state.item(8)
        e3 = state.item(9)
        p = state.item(10)
        q = state.item(11)
        r = state.item(12)
        #   extract forces/moments
        fx = forces_moments.item(0)
        fy = forces_moments.item(1)
        fz = forces_moments.item(2)
        l = forces_moments.item(3)
        m = forces_moments.item(4)
        n = forces_moments.item(5)

        # position kinematics
        # pos_dot = 
        north_dot = (e1**2 + e0**2 - e2**2 - e3**2) * u + (2 * (e1 * e2 - e3 * e0)) * v + (2 * (e1 * e3 + e2 * e0)) * w
        east_dot = (2 * (e1 * e2 + e3 * e0)) * u + (e2**2 + e0**2 - e1**2 - e3**2) * v + (2 * (e2 * e3 - e1 * e0)) * w
        down_dot = (2 * (e1 * e3 - e2 * e0)) * u + (2 * (e2 * e3 + e1 * e0)) * v + (e3**2 + e0**2 - e1**2 - e2**2) * w

        # position dynamics
        u_dot = r * v - q * w + (1 / MAV.mass) * fx
        v_dot = p * w - r * u + (1 / MAV.mass) * fy
        w_dot = q * u - p * v + (1 / MAV.mass) * fz

        # rotational kinematics
        e0_dot = 1/2 * (0 * e0 + -p * e1 + -q * e2 + -r * e3)
        e1_dot = 1/2 * (p * e0 + 0 * e1 + r * e2 + -q * e3)
        e2_dot = 1/2 * (q * e0 + -r * e1 + 0 * e2 + p * e3)
        e3_dot = 1/2 * (r * e0 + q * e1 + -p * e2 + 0 * e3)

        # rotational dynamics
        p_dot = MAV.gamma1 * p * q - MAV.gamma2 * q * r + MAV.gamma3 * l + MAV.gamma4 * n
        q_dot = MAV.gamma5 * p * r - MAV.gamma6 * (p**2 - r**2) + (1 / MAV.Jy) * m
        r_dot = MAV.gamma7 * p * q - MAV.gamma1 * q * r + MAV.gamma4 * l + MAV.gamma8 * n

        # collect the derivative of the states
        x_dot = np.array([[north_dot, east_dot, down_dot, u_dot, v_dot, w_dot,
                           e0_dot, e1_dot, e2_dot, e3_dot, p_dot, q_dot, r_dot]]).T
        return x_dot

    def _update_velocity_data(self, wind=np.zeros((6, 1))):
        steady_state = wind[0:3]
        gust = wind[3:6]

        # convert wind vector from world to body frame and add gust
        wind_body_frame = Quaternion2Rotation(self._state[6:10]) @ steady_state + gust

        # velocity vector relative to the airmass
        v_air = self._state[3:6] - wind_body_frame
        ur = v_air[0]
        vr = v_air[1]
        wr = v_air[2]

        # compute airspeed
        self._Va = (ur**2 + vr**2 + wr**2) ** (1/2)

        # compute angle of attack
        if ur == 0:
            self._alpha = np.pi/2
        else:
            self._alpha = np.arctan(wr/ur)

        # compute sideslip angle
        if self._Va == 0:
            self._beta = 0
        else:
            self._beta = np.arcsin(vr/((ur**2 + vr**2 + wr**2) ** (1/2)))

    def _forces_moments(self, delta):
        """
        return the forces on the UAV based on the state, wind, and control surfaces
        :param delta: np.matrix(delta_a, delta_e, delta_r, delta_t)
        :return: Forces and Moments on the UAV np.matrix(Fx, Fy, Fz, Ml, Mn, Mm)
        """

        phi, theta, psi = Quaternion2Euler(self._state[6:10])
        p = self._state.item(10)
        q = self._state.item(11)
        r = self._state.item(12)

        # compute gravitational forces
        fg = Quaternion2Rotation(self._state[6:10]).T @ np.array([[0], [0], [MAV.mass * MAV.gravity]])
        e_neg_M = math.exp(-MAV.M * (self._alpha - MAV.alpha0))
        e_pos_M = math.exp(MAV.M * (self._alpha + MAV.alpha0))
        sigma = (1 + e_neg_M + e_pos_M) / ((1 + e_neg_M) * (1 + e_pos_M))
        # compute Lift and Drag coefficients
        # non-linear model
        # CL = (np.pi * MAV.AR) / (1 + ((1 + (MAV.AR / 2)**2) ** (1/2))) # THIS IS C L SUB ALPHA (Just before 4.11)
        # CD = MAV.C_D_p + ((MAV.C_L_0 + CL * self._alpha)**2 / (np.pi * MAV.e * MAV.AR)) # THIS IS C D OF ALPHA (4.11) (Non-linear)

        CL = (1 - sigma) * (MAV.C_L_0 + MAV.C_L_alpha * self._alpha) + sigma * (
                    2 * np.sign(self._alpha) * (np.sin(self._alpha) ** 2) * np.cos(self._alpha))
        CD = MAV.C_D_p + ((MAV.C_L_0 + MAV.C_L_alpha * self._alpha) ** 2) / (np.pi * MAV.e * MAV.AR)

        # compute Lift and Drag Forces
        # p 45
        # F_lift = 1/2 * MAV.rho * self._Va**2 * MAV.S_wing * (MAV.C_L_0 +
        #         MAV.C_L_alpha * self._alpha + MAV.C_L_q * ((MAV.c * q) / (2 * self._Va)) +
        #         MAV.C_L_delta_e * delta.elevator)

        # F_drag = 1/2 * MAV.rho * self._Va**2 * MAV.S_wing * (MAV.C_D_0 +
        #         MAV.C_D_alpha * self._alpha + MAV.C_D_q * ((MAV.c * q) / (2 * self._Va)) +
        #         MAV.C_D_delta_e * delta.elevator)

        F_lift = 1/2 * MAV.rho * self._Va**2 * MAV.S_wing * (CL + MAV.C_L_q * MAV.c / (2 * self._Va) * q + MAV.C_L_delta_e * delta.elevator)
        F_drag = 1/2 * MAV.rho * self._Va**2 * MAV.S_wing * (CD + MAV.C_D_q * MAV.c / (2 * self._Va) * q + MAV.C_D_delta_e * delta.elevator)
        # compute propeller thrust and torque
        thrust_prop, torque_prop = self._motor_thrust_torque(self._alpha, delta.throttle)
        # compute propeller thrust and torque
        # thrust_prop, torque_prop = self._motor_thrust_torque(self._Va, delta.throttle)

        # compute longitudinal forces in body frame
        # p 49

        # You can screw yourself I typed all of this out nicely, I refuse to use F_list and F_drag and matrix
        # multiply this out :) -BA
        # CL_of_alpha = MAV.C_L_0 + MAV.C_L_alpha*self._alpha
        # CD_of_alpha = MAV.C_D_0 + MAV.C_D_alpha*self._alpha
        q_dynamic = 0.5 * MAV.rho * (self._Va ** 2) * MAV.S_wing
        fx = -F_drag * np.cos(self._alpha) + F_lift * np.sin(self._alpha) + fg.item(0) + thrust_prop
        fz = -F_drag * np.sin(self._alpha) - F_lift * np.cos(self._alpha) + fg.item(2)
        # compute lateral forces in body frame
        b_2Va = MAV.b / (2 * self._Va)
        beta = self.true_state.beta
        fy = 1/2 * MAV.rho * self._Va**2 * MAV.S_wing * (
                    MAV.C_Y_0 + MAV.C_Y_beta * beta + MAV.C_Y_p * b_2Va * p + MAV.C_Y_r * b_2Va * r + MAV.C_Y_delta_a * delta.aileron + MAV.C_Y_delta_r * delta.rudder) + fg.item(
            1)
        # compute logitudinal torque in body frame
        My = 1/2 * MAV.rho * self._Va**2 * MAV.S_wing * MAV.c * (
                    MAV.C_m_0 + MAV.C_m_alpha * self._alpha + MAV.C_m_q * MAV.c / (2 * self._Va) * q + MAV.C_m_delta_e * delta.elevator)

        # compute lateral torques in body frame
        Mx = 1/2 * MAV.rho * self._Va**2 * MAV.S_wing * MAV.b * (
                    MAV.C_ell_0 + MAV.C_ell_beta * beta + MAV.C_ell_p * b_2Va * p + MAV.C_ell_r * b_2Va * r + MAV.C_ell_delta_a * delta.aileron + MAV.C_ell_delta_r * delta.rudder) + torque_prop
        Mz = 1/2 * MAV.rho * self._Va**2 * MAV.S_wing * MAV.b * (
                    MAV.C_n_0 + MAV.C_n_beta * beta + MAV.C_n_p * b_2Va * p + MAV.C_n_r * b_2Va * r + MAV.C_n_delta_a * delta.aileron + MAV.C_n_delta_r * delta.rudder)
        
        # fx = 1/2 * MAV.rho * self._Va**2 * MAV.S_wing * (
        #         (-CD_of_alpha * np.cos(self._alpha) + CL_of_alpha * np.sin(self._alpha)) +
        #         (-MAV.C_D_q * np.cos(self._alpha) + MAV.C_L_q * np.sin(self._alpha)) * ((MAV.c * q) / (2 * self._Va)) +
        #         (-MAV.C_D_delta_e * np.cos(self._alpha) + MAV.C_L_delta_e * np.sin(self._alpha)) * delta.elevator)

        # fz = 1/2 * MAV.rho * self._Va**2 * MAV.S_wing * (
        #         (-CD_of_alpha * np.sin(self._alpha) - CL_of_alpha * np.cos(self._alpha)) +
        #         (-MAV.C_D_q * np.sin(self._alpha) - MAV.C_L_q * np.cos(self._alpha)) * ((MAV.c * q) / (2 * self._Va)) +
        #         (-MAV.C_D_delta_e * np.sin(self._alpha) - MAV.C_L_delta_e * np.cos(self._alpha)) * delta.elevator)

        # compute lateral forces in body frame
        # p 50
        # fy = 1/2 * MAV.rho * self._Va**2 * MAV.S_wing * (
        #         (MAV.C_Y_0 + MAV.C_Y_beta * self._beta + MAV.C_Y_p * (MAV.b * p)/(2 * self._Va) +
        #         MAV.C_Y_r * (MAV.b * r)/(2 * self._Va) + MAV.C_Y_delta_a * delta.aileron + MAV.C_Y_delta_r * delta.rudder))

        # # compute longitudinal torque in body frame
        # # p 45
        # My = (MAV.rho / 2) * self._Va**2 * MAV.S_wing * MAV.c * (MAV.C_m_0 + MAV.C_m_alpha * self._alpha +
        #         MAV.C_m_q * (MAV.c * q) / (2 * self._Va) + MAV.C_m_delta_e * delta.elevator)
                # Pitching moment

        # compute lateral torques in body frame
        # p 50-51
        # 4.15
        # Mx = (MAV.rho/2) * self._Va**2 * MAV.S_wing * MAV.b * (MAV.C_ell_0 + MAV.C_ell_beta * self._beta +
        #         MAV.C_ell_p * (MAV.b * p)/(2 * self._Va) + MAV.C_ell_r * (MAV.b * r)/(2 * self._Va) +
        #         MAV.C_ell_delta_a * delta.aileron + MAV.C_ell_delta_r * delta.rudder)
        #         # Roll moment
        
        # # 4.16
        # Mz = (MAV.rho/2) * self._Va**2 * MAV.S_wing * MAV.b * (MAV.C_n_0 + MAV.C_n_beta * self._beta +
        #         MAV.C_n_p * (MAV.b * p)/(2 * self._Va) + MAV.C_n_r * (MAV.b * r)/(2 * self._Va) +
        #         MAV.C_n_delta_a * delta.aileron + MAV.C_n_delta_r * delta.rudder)
                # Yaw moment

        self._forces[0] = fx
        self._forces[1] = fy
        self._forces[2] = fz

        return np.array([[fx, fy, fz, Mx, My, Mz]]).T

    def _motor_thrust_torque(self, Va, delta_t):
        # compute thrust and torque due to propeller  (See addendum by McLain)
        # map delta_t throttle command(0 to 1) into motor input voltage
        V_Max = 3.7 * MAV.ncells
        V_in = V_Max * delta_t

        # Angular speed of propeller
        # slide 34
        a = (MAV.rho * MAV.D_prop**5 * MAV.C_Q0) / ((2*np.pi)**2)
        b = ((MAV.rho * MAV.D_prop**4 * MAV.C_Q1 * self._Va) / (2*np.pi)) + (MAV.KQ * MAV.KV) / MAV.R_motor
        c = MAV.rho * MAV.D_prop**3 * MAV.C_Q2 * self._Va**2 - MAV.KQ * V_in / MAV.R_motor + MAV.KQ * MAV.i0
        Omega_p = (-b + (b**2 - 4*a*c) ** (0.5)) / (2.*a)

        J_op = 2 * np.pi * self._Va / (Omega_p * MAV.D_prop)

        C_T = MAV.C_T2 * J_op **2 + MAV.C_T1 * J_op + MAV.C_T0
        C_Q = MAV.C_Q2 * J_op **2 + MAV.C_Q1 * J_op + MAV.C_Q0
        n = Omega_p / (2*np.pi)

        # thrust and torque due to propeller
        # slide 31

        #get from slides chapter 4 
        thrust_prop = MAV.rho * n**2 * MAV.D_prop**4 * C_T
        torque_prop = -MAV.rho * n**2 * MAV.D_prop**5 * C_Q 
        return thrust_prop, torque_prop

    def _update_true_state(self):
        # update the class structure for the true state:
        #   [pn, pe, h, Va, alpha, beta, phi, theta, chi, p, q, r, Vg, wn, we, psi, gyro_bx, gyro_by, gyro_bz]
        phi, theta, psi = Quaternion2Euler(self._state[6:10])
        pdot = Quaternion2Rotation(self._state[6:10]) @ self._state[3:6]
        self.true_state.north = self._state.item(0)
        self.true_state.east = self._state.item(1)
        self.true_state.altitude = -self._state.item(2)
        self.true_state.Va = self._Va
        self.true_state.alpha = self._alpha
        self.true_state.beta = self._beta
        self.true_state.phi = phi
        self.true_state.theta = theta
        self.true_state.psi = psi
        self.true_state.Vg = np.linalg.norm(pdot)
        self.true_state.gamma = np.arcsin(pdot.item(2) / self.true_state.Vg)
        self.true_state.chi = np.arctan2(pdot.item(1), pdot.item(0))
        self.true_state.p = self._state.item(10)
        self.true_state.q = self._state.item(11)
        self.true_state.r = self._state.item(12)
        self.true_state.wn = self._wind.item(0)
        self.true_state.we = self._wind.item(1)
