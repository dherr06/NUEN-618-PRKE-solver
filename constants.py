import numpy as np

class constants:
    
    def P():
        P = 59506.0 #initial power in single pin
        return P
    def gamma():
        gamma = np.array([0.1])
        return gamma
    def beta():
        beta = np.array([600e-5])
        return beta
    def zeta():
        zeta = np.array([constants.beta()[0]*constants.P()/(constants.gamma()[0])])
        return zeta
    def eta():
        eta = 1e-5
        return eta
    def source():
        source = 0.0
        return source
    def alpha_f():
        alpha_f = -3.5e-5
        return alpha_f
    def alpha_m():
        alpha_m = -55e-5
        return alpha_m
    def R_f():
        R_f = 0.0041
        return R_f
    def R_g():
        R_g = 0.00411
        return R_g
    def R_c():
        R_c = 0.00475
        return R_c
    def fuel_height():
        fuel_height = 4
        return fuel_height
    def pitch():
        pitch = 0.0126
        return pitch
    def gap_cond():
        gap_cond = 10**4
        return gap_cond
    def clad_cond():
        clad_cond = 17
        return clad_cond
    def fluid_axial_velocity():
        fluid_axial_velocity = 5
        return fluid_axial_velocity
    def fluid_cond():
        fluid_cond = 0.54
        return fluid_cond
    def fluid_visc():
        fluid_visc = 90*1e-6
        return fluid_visc
    def T_inlet():
        T_inlet = 290
        return T_inlet
    def T_cool_SS():
        T_cool_SS = 306.8830857
        return T_cool_SS
    def T_fuel_SS():
        T_fuel_SS = 742.3225084
        return T_fuel_SS