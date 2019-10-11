import numpy as np
from scipy.optimize import fsolve
from properties import ragusa_code as RC
from constants import constants
    
class PRKE:
    
    def crank_nic(P,zeta,rho_0,rho_1,t,eta,gamma,beta,source):
        '''
        This functions solves the PRKE's using crank nicholson
        '''
        #initialize containers
        A = np.zeros((np.size(beta)+1,np.size(beta)+1))
        b = np.zeros(np.size(beta)+1)
        #compute the total delayed neutron fraction
        beta_total = np.sum(beta)
        #fill the matrix and the b vector
        A[0][0] = ((t/2)*((rho_1 - beta_total)/eta)) - 1.0
        b[0] = -P - ((t/2) * ((rho_0 - beta_total)/eta) * P) - (t*source/eta)
        for i in range(1,np.size(beta)+1):
            A[0][i] = (t/2.0) * (gamma[i-1]/eta)
            b[0] -= ((t/2.0) * (gamma[i-1]/eta) * (zeta[i-1]))
            A[i][0] = t*beta[i-1]/2.0
            A[i][i] = -1.0 - (gamma[i-1]*t/2.0)
            #fill b
            b[i] = -zeta[i-1] + ((t/2)*gamma[i-1]*zeta[i-1]) - ((t/2)*beta[i-1]*P)
        #solve the linear system
        soln = np.linalg.solve(A,b)
        return soln
    
class non_linear:
    
    def Newton(x,func,epsilon):
        '''
        Performs a non linear solve using a finite difference Jacobian
        '''
        J = np.zeros((func.shape[0],x.shape[0]))
        for i in range(x.shape[0]):
            e = np.zeros(x.shape[0])
            e[i] = 1
            J[:,i] = (func(x + e*epsilon) - func(x))/epsilon
        R = np.linalg.solve(J,func(x))
        x_new = x - R
        
class G:
    
    def __init__(self,rho_0,t,P,zeta,T_cool,T_fuel,staggered=True):
        self.rho_0 = rho_0
        self.t = t
        self.T_cool = T_cool
        self.T_fuel = T_fuel
        self.P = P
        self.zeta = zeta
        self.staggered = staggered
    
    def solve(self,rho_1,P_g,T_cool_g,T_fuel_g):
        #initialialize constants
        rho_0 = self.rho_0
        t = self.t
        P = self.P
        zeta = self.zeta
        T_cool = self.T_cool
        T_fuel = self.T_fuel
        staggered = self.staggered
        A_fuel = np.pi*constants.R_f()**2
        A_flow = constants.pitch()**2 - (np.pi*constants.R_c()**2)
        P_wet = 2*np.pi*constants.R_c()
        D_hy = 4*A_flow/P_wet
        #define Re and Pr as functions of temp
        Re = lambda T: RC.rho_mod(T)[0]*constants.fluid_axial_velocity()*D_hy/constants.fluid_visc()
        Pr = lambda T: constants.fluid_visc()*RC.cp_mod(T)[0]/constants.fluid_cond()
        h_conv = lambda T:0.023*(Re(T)**0.8)*(Pr(T)**0.4)*(constants.fluid_cond()/D_hy)
        R_th = lambda T_f,T_c: A_fuel/2/np.pi*(1/constants.R_g()/constants.gap_cond() + 1/constants.clad_cond()*np.log(constants.R_c()/constants.R_g()) + 1/constants.R_c()/h_conv(T_c) + 1/2/RC.k_fuel(T_f)[0])
        #if staggered, perform staggered OS on all physics
        if staggered == True:
            P_new,zeta_new = PRKE.crank_nic(P,zeta,rho_0,rho_1,t,constants.eta(),constants.gamma(),constants.beta(),0)
            fuel_func = lambda x: T_fuel - x + (t/2)*(1/RC.rhocp_fuel(x)[0] * (P/constants.fuel_height()/A_fuel - ((x - T_cool)/R_th(x,T_cool))) + 1/RC.rhocp_fuel(x)[0] * (P_new/constants.fuel_height()/A_fuel - ((x - T_cool_g)/R_th(x,T_cool_g))))
            T_fuel_new = fsolve(fuel_func,T_fuel)
            coolant_func = lambda x: T_cool - x + (t/2)*((1/RC.rhocp_mod(x)[0]/A_flow * A_fuel*((T_fuel - x)/R_th(T_fuel,x)) - (constants.fluid_axial_velocity()*2/constants.fuel_height()*(x - constants.T_inlet()))) + (1/RC.rhocp_mod(x)[0]/A_flow * A_fuel*((T_fuel_new - x)/R_th(T_fuel_new,x)) - (constants.fluid_axial_velocity()*2/constants.fuel_height()*(x - constants.T_inlet()))))
            T_cool_new = fsolve(coolant_func,T_cool)
        #perform simultaenous solves on all physics
        else:
            P_new,zeta_new = PRKE.crank_nic(P,zeta,rho_0,rho_1,t,constants.eta(),constants.gamma(),constants.beta(),0)
            def funcs(T_f,T_c):
                fuel_func = T_fuel - T_f + (t/2)*(1/RC.rhocp_fuel(T_f)[0] * (P/constants.fuel_height()/A_fuel - ((T_f - T_cool)/R_th(T_f,T_cool))) + 1/RC.rhocp_fuel(T_f)[0] * (P_g/constants.fuel_height()/A_fuel - ((T_f - T_cool_g)/R_th(T_f,T_cool_g))))
                coolant_func = T_cool - T_c + (t/2)*((1/RC.rhocp_mod(T_c)[0]/A_flow * A_fuel*((T_fuel - T_c)/R_th(T_fuel,T_c)) - (constants.fluid_axial_velocity()*2/constants.fuel_height()*(T_c - constants.T_inlet()))) + (1/RC.rhocp_mod(T_c)[0]/A_flow * A_fuel*((T_fuel_g - x)/R_th(T_fuel_g,T_c)) - (constants.fluid_axial_velocity()*2/constants.fuel_height()*(x - constants.T_inlet()))))
                return fuel_func,coolant_func
            if mode == 'fsolve':
                T_fuel_new, T_cool_new = fsolve(funcs,(T_fuel,T_cool))
            elif mode == 'newton':
                T_fuel_new = non_linear.Newton(x,funcs,1e-12)
        #return the end time values
        return P_new, np.array([zeta_new]), T_cool_new, T_fuel_new
