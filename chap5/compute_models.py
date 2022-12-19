"""
compute_ss_model
    - Chapter 5 assignment for Beard & McLain, PUP, 2012
    - Update history:  
        2/4/2019 - RWB
"""
import sys
sys.path.append('..')
import numpy as np
from scipy.optimize import minimize
from tools.rotations import Euler2Quaternion, Quaternion2Euler
import parameters.aerosonde_parameters as MAV
from parameters.simulation_parameters import ts_simulation as Ts
from message_types.msg_delta import MsgDelta
from control.matlab import *

def compute_model(mav, trim_state, trim_input):
    A_lon, B_lon, A_lat, B_lat = compute_ss_model(mav, trim_state, trim_input)
    Va_trim, alpha_trim, theta_trim, a_phi1, a_phi2, a_theta1, a_theta2, a_theta3, \
    a_V1, a_V2, a_V3 = compute_tf_model(mav, trim_state, trim_input)

    # write transfer function gains to file
    file = open('model_coef.py', 'w')
    file.write('import numpy as np\n')
    file.write('x_trim = np.array([[%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f]]).T\n' %
               (trim_state.item(0), trim_state.item(1), trim_state.item(2), trim_state.item(3),
                trim_state.item(4), trim_state.item(5), trim_state.item(6), trim_state.item(7),
                trim_state.item(8), trim_state.item(9), trim_state.item(10), trim_state.item(11),
                trim_state.item(12)))
    file.write('u_trim = np.array([[%f, %f, %f, %f]]).T\n' %
               (trim_input.elevator, trim_input.aileron, trim_input.rudder, trim_input.throttle))
    file.write('Va_trim = %f\n' % Va_trim)
    file.write('alpha_trim = %f\n' % alpha_trim)
    file.write('theta_trim = %f\n' % theta_trim)
    file.write('a_phi1 = %f\n' % a_phi1)
    file.write('a_phi2 = %f\n' % a_phi2)
    file.write('a_theta1 = %f\n' % a_theta1)
    file.write('a_theta2 = %f\n' % a_theta2)
    file.write('a_theta3 = %f\n' % a_theta3)
    file.write('a_V1 = %f\n' % a_V1)
    file.write('a_V2 = %f\n' % a_V2)
    file.write('a_V3 = %f\n' % a_V3)
    file.write('A_lon = np.array([\n    [%f, %f, %f, %f, %f],\n    '
               '[%f, %f, %f, %f, %f],\n    '
               '[%f, %f, %f, %f, %f],\n    '
               '[%f, %f, %f, %f, %f],\n    '
               '[%f, %f, %f, %f, %f]])\n' %
    (A_lon[0][0], A_lon[0][1], A_lon[0][2], A_lon[0][3], A_lon[0][4],
     A_lon[1][0], A_lon[1][1], A_lon[1][2], A_lon[1][3], A_lon[1][4],
     A_lon[2][0], A_lon[2][1], A_lon[2][2], A_lon[2][3], A_lon[2][4],
     A_lon[3][0], A_lon[3][1], A_lon[3][2], A_lon[3][3], A_lon[3][4],
     A_lon[4][0], A_lon[4][1], A_lon[4][2], A_lon[4][3], A_lon[4][4]))
    file.write('B_lon = np.array([\n    [%f, %f],\n    '
               '[%f, %f],\n    '
               '[%f, %f],\n    '
               '[%f, %f],\n    '
               '[%f, %f]])\n' %
    (B_lon[0][0], B_lon[0][1],
     B_lon[1][0], B_lon[1][1],
     B_lon[2][0], B_lon[2][1],
     B_lon[3][0], B_lon[3][1],
     B_lon[4][0], B_lon[4][1],))
    file.write('A_lat = np.array([\n    [%f, %f, %f, %f, %f],\n    '
               '[%f, %f, %f, %f, %f],\n    '
               '[%f, %f, %f, %f, %f],\n    '
               '[%f, %f, %f, %f, %f],\n    '
               '[%f, %f, %f, %f, %f]])\n' %
    (A_lat[0][0], A_lat[0][1], A_lat[0][2], A_lat[0][3], A_lat[0][4],
     A_lat[1][0], A_lat[1][1], A_lat[1][2], A_lat[1][3], A_lat[1][4],
     A_lat[2][0], A_lat[2][1], A_lat[2][2], A_lat[2][3], A_lat[2][4],
     A_lat[3][0], A_lat[3][1], A_lat[3][2], A_lat[3][3], A_lat[3][4],
     A_lat[4][0], A_lat[4][1], A_lat[4][2], A_lat[4][3], A_lat[4][4]))
    file.write('B_lat = np.array([\n    [%f, %f],\n    '
               '[%f, %f],\n    '
               '[%f, %f],\n    '
               '[%f, %f],\n    '
               '[%f, %f]])\n' %
    (B_lat[0][0], B_lat[0][1],
     B_lat[1][0], B_lat[1][1],
     B_lat[2][0], B_lat[2][1],
     B_lat[3][0], B_lat[3][1],
     B_lat[4][0], B_lat[4][1],))
    file.write('Ts = %f\n' % Ts)
    file.close()


def compute_tf_model(mav, trim_state, trim_input):
    # trim values
    mav._state = trim_state
    mav._update_velocity_data()
    Va_trim = mav._Va
    alpha_trim = mav._alpha
    phi, theta_trim, psi = Quaternion2Euler(trim_state[6:10])
    delta_trim = mav._delta

    # define transfer function constants
    a_phi1 = -1/2 * MAV.rho * (Va_trim ** 2) *  MAV.S_wing * MAV.b * MAV.C_p_p * MAV.b / (2 * Va_trim)
    a_phi2 = 1/2 * MAV.rho * (Va_trim ** 2) *  MAV.S_wing * MAV.b * MAV.C_p_delta_a
    a_theta1 = -((MAV.rho * Va_trim ** 2 * MAV.c * MAV.S_wing) / (2 * MAV.Jy)) * MAV.C_m_q
    a_theta2 = -((MAV.rho * Va_trim ** 2 * MAV.c * MAV.S_wing) / (2 * MAV.Jy)) * MAV.C_m_alpha
    a_theta3 = ((MAV.rho * Va_trim ** 2 * MAV.c * MAV.S_wing) / (2 * MAV.Jy)) * MAV.C_m_delta_e

    # Compute transfer function coefficients using new propulsion model
    # CHAPTER 5: SLIDE 19: LEC 7-8 LINEAR MODELS
    
    # REVIEW AV1
    # PAGE 94 TB?
    a_V1 = ((MAV.rho * Va_trim * MAV.S_wing) / MAV.mass) * \
           (MAV.C_D_0 + MAV.C_D_alpha * alpha_trim + MAV.C_D_delta_e * delta_trim) - \
            (MAV.rho * MAV.S_prop * MAV.C_prop * Va_trim / MAV.mass)
           #dT_dVa(delta_trim, Va_trim) / MAV.mass # AT : Not sure where I got this from (slides?) but the dt_dtva stuff is incorrect
    #a_V2 = dT_ddelta_t(mav, delta_trim, Va_trim) / MAV.mass
    a_V2 = (MAV.rho * MAV.S_prop * MAV.C_prop * MAV.k_motor ** 2 ) / MAV.mass
    a_V3 = MAV.gravity * np.cos(theta_trim - alpha_trim)

    return Va_trim, alpha_trim, theta_trim, a_phi1, a_phi2, a_theta1, a_theta2, a_theta3, a_V1, a_V2, a_V3


def compute_ss_model(mav, trim_state, trim_input):
    # in chpt 5
    mav._state = trim_state
    mav._update_velocity_data()
    x_euler = euler_state(trim_state)
    
    pn = x_euler.item(0)
    pe = x_euler.item(1)
    pd = x_euler.item(2)
    u = x_euler.item(3)
    v = x_euler.item(4)
    w = x_euler.item(5)
    phi = x_euler.item(6)
    theta = x_euler.item(7)
    psi = x_euler.item(8)
    p = x_euler.item(9)
    q = x_euler.item(10)
    r = x_euler.item(11)
    delta_a = trim_input.item(0)
    delta_e = trim_input.item(1)
    delta_r = trim_input.item(2)
    delta_t = trim_input.item(3)

    beta_star = mav._beta
    Va_star = mav._Va
    alpha_trim = mav._alpha
    u_star = u
    v_star = v
    w_star = w
    phi_star = phi
    theta_star = theta
    p_star = p
    q_star = q
    r_star = r

    delta_a_star = delta_a
    delta_r_star = delta_r
    

    CD = MAV.C_D_p + ((MAV.C_L_0 + alpha_trim * MAV.C_L_alpha)**2)/(np.pi * MAV.e * MAV.AR)
    Sigma = (1.0 + np.exp(-MAV.M * (alpha_trim - MAV.alpha0))+np.exp(MAV.M*(alpha_trim + MAV.alpha0))) / (1.0 + np.exp(-MAV.M * (alpha_trim - MAV.alpha0))) * (1 + np.exp(MAV.M * (alpha_trim + MAV.alpha0)))
    CL = (1-Sigma)*(MAV.C_L_0 + MAV.C_L_alpha*alpha_trim)+Sigma*(2.0*np.sign(alpha_trim)*np.sin(alpha_trim)**2*np.cos(alpha_trim))
    
    Cx = -CD * np.cos(alpha_trim) + CL * np.sin(alpha_trim)
    
    Xu = u_star*MAV.rho*MAV.S_wing*(Cx)
    Xw = 0.0
    Xq = 0.0
    X_delta_e = 0.0
    X_delta_t = 0.0
    Zu = 0.0
    Zw = 0.0
    Zq = 0.0
    Z_delta_e = 0.0
    Mu = 0.0
    Mw = 0.0
    Mq = 0.0
    M_delta_e = 0.0
    
    #A = 
    #B = 
    # extract longitudinal states (u, w, q, theta, pd) and change pd to h
    A_lon = np.array([[Xu, Xw * Va_star * np.cos(alpha_trim), Xq, -MAV.gravity * np.cos(theta_star),0.0],
                      [Zu / (Va_star * np.cos(alpha_trim)),Zw,Zq/(Va_star*np.cos(alpha_trim)), -MAV.gravity * np.sin(theta_star) / (Va_star * np.cos(alpha_trim)),0.0],
                      [Mu, Mw * Va_star, np.cos(alpha_trim), Mq,0.0,0.0],
                      [0.0,0.0,1.0,0.0,0.0],
                      [np.sin(theta_star), -Va_star * np.cos(theta_star) * np.cos(alpha_trim),0.0,u_star * np.cos(theta_star)+w_star * np.sin(theta_star),0.0]])
    B_lon = np.array([[X_delta_e, X_delta_t],
                      [Z_delta_e / (Va_star * np.cos(alpha_trim)) ,0.0],
                      [M_delta_e, 0.0],
                      [0.0,0.0],
                      [0.0,0.0]])

    Yv = MAV.rho * MAV.S_wing * MAV.b*v_star * (MAV.C_Y_p*p_star+MAV.C_Y_r * r_star)/(4.*MAV.mass * Va_star) \
        + MAV.rho * MAV.S_wing * v_star * (MAV.C_Y_0+MAV.C_Y_beta*beta_star+MAV.C_Y_delta_a*delta_a_star+MAV.C_Y_delta_r*delta_r_star)/MAV.mass \
        + MAV.rho * MAV.S_wing * MAV.C_Y_beta * np.sqrt(u_star**2+w_star**2)/(2.*MAV.mass)
    
    Yp = w_star + MAV.rho * Va_star*MAV.S_wing*MAV.b*MAV.C_Y_p/(4.*MAV.mass)
    
    Yr = -u_star + MAV.rho*Va_star*MAV.S_wing*MAV.b*MAV.C_Y_r/(4.*MAV.mass)
    
    Y_delta_a = MAV.rho*Va_star**2*MAV.S_wing*MAV.C_Y_delta_a/(2.*MAV.mass)
    
    Y_delta_r = MAV.rho*Va_star**2*MAV.S_wing*MAV.C_Y_delta_r/(2.*MAV.mass)
    
    Lv = MAV.rho*MAV.S_wing*MAV.b**2*v_star*(MAV.C_p_p*p_star+MAV.C_p_r*r_star)/(4.*Va_star) \
        + MAV.rho*MAV.S_wing*MAV.b*v_star*(MAV.C_p_0+MAV.C_p_beta*beta_star+MAV.C_p_delta_a*delta_a_star+MAV.C_p_delta_r*delta_r_star) \
        + MAV.rho*MAV.S_wing*MAV.b*MAV.C_p_beta*np.sqrt(u_star**2+w_star**2)/(2.)
    
    Lp = MAV.gamma1*q_star + MAV.rho*Va_star*MAV.S_wing*MAV.b**2*MAV.C_p_p/4.
    
    Lr = -MAV.gamma2*q_star + MAV.rho*Va_star*MAV.S_wing*MAV.b**2*MAV.C_p_r/4.
    
    L_delta_a = MAV.rho*Va_star**2*MAV.S_wing*MAV.b*MAV.C_p_delta_a/2.
    
    L_delta_r = MAV.rho*Va_star**2*MAV.S_wing*MAV.b*MAV.C_p_delta_r/2.
    
    Nv = MAV.rho*MAV.S_wing*MAV.b**2*v_star*(MAV.C_r_p*p_star+MAV.C_r_r*r_star)/(4.*Va_star) \
        + MAV.rho*MAV.S_wing*MAV.b*v_star*(MAV.C_r_0+MAV.C_r_beta*beta_star+MAV.C_r_delta_a*delta_a_star+MAV.C_r_delta_r*delta_r_star) \
        + MAV.rho*MAV.S_wing*MAV.b*MAV.C_r_beta*np.sqrt(u_star**2+w_star**2)/(2.)
    
    Np = MAV.gamma7*q_star + MAV.rho*Va_star*MAV.S_wing*MAV.b**2*MAV.C_r_p/4.
    
    Nr = -MAV.gamma1*q_star + MAV.rho*Va_star*MAV.S_wing*MAV.b**2*MAV.C_r_r/4.
    
    N_delta_a = MAV.rho*Va_star**2*MAV.S_wing*MAV.b*MAV.C_r_delta_a/2.
    
    N_delta_r = MAV.rho*Va_star**2*MAV.S_wing*MAV.b*MAV.C_r_delta_r/2.

    # extract lateral states (v, p, r, phi, psi)
    A_lat = np.array([[Yv, Yp/(Va_star*np.cos(beta_star)), Yr/(Va_star*np.cos(beta_star)), MAV.gravity*np.cos(theta_star)*np.cos(phi_star)/(Va_star*np.cos(beta_star)),0.0],
                      [Lv*Va_star*np.cos(beta_star), Lp, Lr, 0.0, 0.0],
                      [Nv*Va_star*np.cos(beta_star),Np,Nr,0.0,0.0],
                      [0.0, 1.0, np.cos(phi_star)*np.tan(theta_star),q_star*np.cos(phi_star)*np.tan(theta_star)-r_star*np.sin(phi_star)*np.tan(theta_star),0.0],
                      [0.0, 0.0, np.cos(phi_star)/np.cos(theta_star),p_star*np.cos(phi_star)/np.cos(theta_star)-r_star*np.sin(phi_star)/np.cos(theta_star),0.0]])
    B_lat = np.array([[Y_delta_a/(Va_star*np.cos(beta_star)),Y_delta_r/(Va_star*np.cos(beta_star))],
                      [L_delta_a,L_delta_r],
                      [N_delta_a,N_delta_r],
                      [0.0, 0.0],
                      [0.0, 0.0]])
    return A_lon, B_lon, A_lat, B_lat

def euler_state(x_quat):
    # convert state x with attitude represented by quaternion
    # to x_euler with attitude represented by Euler angles

    # THIS MIGHT BE VERY WRONG, BUT I THINK IT'S THE RIGHT IDEA

    # BA - I think you're good here. Did a bit of redefining for clarity.

    e = np.array([x_quat.item(6), x_quat.item(7), x_quat.item(8), x_quat.item(9)])
    angles = np.array(Quaternion2Euler(e))
    #angles = Quaternion2Euler(e)
    # print("angles:",  angles)
    # print("angles[0]:",  angles[0])

    x_euler = np.array([[x_quat.item(0)],
                        [x_quat.item(1)],
                        [x_quat.item(2)],
                        [x_quat.item(3)],
                        [x_quat.item(4)],
                        [x_quat.item(5)],
                        [angles[0]], # AT changed these to angles[0], instead of phi, theta, psi (undefined vars)
                        [angles[1]],
                        [angles[2]],
                        [x_quat.item(10)],
                        [x_quat.item(11)],
                        [x_quat.item(12)]
                        ])
    return x_euler

def quaternion_state(x_euler):
    # convert state x_euler with attitude represented by Euler angles
    # to x_quat with attitude represented by quaternions

    e = Euler2Quaternion(x_euler.item(6),x_euler.item(7), x_euler.item(8))
    x_quat = np.array([[x_euler.item(0)],
                        [x_euler.item(1)],
                        [x_euler.item(2)],
                        [x_euler.item(3)],
                        [x_euler.item(4)],
                        [x_euler.item(5)],
                        [e.item(0)],
                        [e.item(1)],
                        [e.item(2)],
                        [e.item(3)],
                        [x_euler.item(10)],
                        [x_euler.item(11)],
                        [x_euler.item(12)]
                        ])

    return x_quat

def f_euler(mav, x_euler, delta):
    # return 12x1 dynamics (as if state were Euler state)
    # compute f at euler_state

    xq = quaternion_state(x_euler)
    mav._state = xq
    mav._update_velocity_data()
    fm = mav._forces_moments(delta)
    fq = mav._derivatives(xq, fm)
    fq[2, 0] = -1 * fq[2, 0] # BA - because of that weird p_down dot vs h dot thing

    q0 = xq[6:10]
    Aq = np.zeros((3,4))
    for i in range(4):
        eps = 0.001
        q_eps = np.zeros_like(q0)
        q_eps[i, 0] = eps
        q_eps = q_eps + q0
        Aq[:, i] = ((np.array([Quaternion2Euler(q_eps)]).T - np.array([Quaternion2Euler(q0)]).T) / eps)[:, 0]
    f_euler_ = np.stack(([fq[0:6], np.matmul(Aq, fq[6:10]), fq[10:13]]))
    return f_euler_

# def df_dx(mav, x_euler, delta):
#     # take partial of f_euler with respect to x_euler
#     # PAGE 80, 84?
#     eps =
#     return A


# def df_du(mav, x_euler, delta):
#     # take partial of f_euler with respect to input
#     # PAGE 80, 84?
#     B = 
#     return B

# 
# def dT_dVa(mav, Va, delta_t):
#     # returns the derivative of motor thrust with respect to Va
#     eps = 
#     T_eps, Q_eps = #mav._motor_thrust_torque()
#     T, Q = #mav._motor_thrust_torque()
#     return (T_eps - T) / eps

# def dT_ddelta_t(mav, Va, delta_t):
#     # returns the derivative of motor thrust with respect to delta_t
#     eps = EPS
#     T_eps, Q_eps = #mav._motor_thrust_torque()
#     T, Q = #mav._motor_thrust_torque()
#     return (T_eps - T) / eps