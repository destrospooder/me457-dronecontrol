import numpy as np
class PIDControl:
    def __init__(self, kp, ki, kd, limit, sigma, Ts, flag=True):
        self.kp = kp # Proportional control gain
        self.ki = ki # Integral control gain
        self.kd = kd # Derivative control gain
        self.limit = limit # The output saturates at this limit
        self.sigma = sigma # dirty derivative bandwidth is 1/sigma 
        self.beta = (2.0*sigma-Ts)/(2.0*sigma+Ts)
        self.Ts = Ts # sample rate 
        self.flag = flag
        # if flag == True, then returns
        # u = kp * error + ki * integral(error) + kd * error_dot.
        # else returns
        # u = kp * error + ki * integral(error) - kd * y_dot. 
        self.y_dot = 0.0 # estimated derivative of y 
        self.y_d1 = 0.0 # Signal y delayed by one sample
        self.error_dot = 0.0 # estimated derivative of error 
        self.error_d1 = 0.0 # Error delayed by one sample
        self.integrator = 0.0 # integrator
    
    def PID(self, y_r, y):
        # Compute the current error
        error = y_r - y
        # integrate error using trapazoidal rule
        self.integrator = self.integrator \
        + (self.Ts/2) * (error + self.error_d1)
        # PID Control
        if self.flag:
            # differentiate error
            self.error_dot = self.beta * self.error_dot \
                + (1-self.beta)/self.Ts * (error - self.error_d1) 
            # PID control
            u_unsat = self.kp*error \
                + self.ki*self.integrator \
                + self.kd*self.error_dot
        else:
            # differentiate y
            self.y_dot = self.beta * self.y_dot \
                + (1-self.beta)/self.Ts * (y - self.y_d1) 
            # PID control
            u_unsat = self.kp*error \
                + self.ki*self.integrator \
                - self.kd*self.y_dot 
        # return saturated control signal
        u_sat = self.saturate(u_unsat) # integrator anti - windup
        if self.ki != 0.0:
            self.integrator = self.integrator \
                + 1.0 / self.ki * (u_sat - u_unsat) 
        # update delayed variables
        self.error_d1 = error 
        self.y_d1 = y
        return u_sat

    def PD(self, y_r, y):
        # Compute the current error
        error = y_r - y 
        # PD Control
        if self.flag:
            # differentiate error
            self.error_dot = self.beta * self.error_dot \
                + (1-self.beta)/self.Ts * (error - self.error_d1)
            # PD control
            u_unsat = self.kp*error \
                + self.kd*self.error_dot
        else:
            # differentiate y
            self.y_dot = self.beta * self.y_dot \
                + (1-self.beta)/self.Ts * (y - self.y_d1)
            # PD control
            u_unsat = self.kp*error \
                - self.kd*self.y_dot
        # return saturated control signal
        u_sat = self.saturate(u_unsat) 
        # update delayed variables
        self.error_d1 = error 
        self.y_d1 = y
        return u_sat
    
    def saturate(self,u):
        if abs(u) > self.limit:
            u = self.limit*np.sign(u)
        return u