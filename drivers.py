import numpy as np
from numerical_toolbox import G
from constants import constants

class drivers:
    
    def driver(t,staggered=True,it_solve=True):
        P = constants.P()
        zeta = constants.zeta()
        power_staggered = np.array([P])
        staggered_it = np.array([])
        rho_0 = 0
        rho_1 = 0
        t = t
        T = 0.05
        num_steps = 15*T/t
        if it_solve == True:
            max_iters = 10
        else:
            max_iters = 1
        T_cool = constants.T_cool_SS()
        T_fuel = constants.T_fuel_SS()
        r_ext = np.array([])
        r = np.array([])
        rho_ext = 0
        withdrawl = False
        #outter most layer, FPI on all the physics
        for step in range(int(num_steps)):
            FPI_object = G(rho_0,t,P,zeta,T_cool,T_fuel,staggered)
            #if the rod has been pulled change rho_ext
            if step*t > 0.05 and step*t <= 0.1:
                rho_ext = (1.45*constants.beta()[0]/0.05)*(step*t - 0.05)
                withdrawl = True
            rho_1 = rho_ext + (constants.alpha_f()*(T_fuel - constants.T_fuel_SS())) + (constants.alpha_m()*(T_cool - constants.T_cool_SS()))
            #do FPI
            for i in range(max_iters):
                num_its = i
                P_new, zeta_new, T_cool_new, T_fuel_new = FPI_object.solve(rho_1,P,T_cool,T_fuel)
                if np.abs(P - P_new) < 1e-8 and np.abs(T_cool - T_cool_new) < 1e-8 and np.abs(T_fuel - T_fuel_new) < 1e-8 and it_solve == True:
                    break
                if it_solve == True:
                    #the change in reactivity resulted in a larger change in power than is tolerable
                    #compute a new end time reactivity and try again
                    #update the latest end time values
                    P = P_new
                    zeta = zeta_new
                    T_cool = T_cool_new
                    T_fuel = T_fuel_new
                    #compute the new reactivity
                    if withdrawl == True:
                        rho_1 = rho_ext + (constants.alpha_f()*(T_fuel_new - constants.T_fuel_SS())) + (constants.alpha_m()*(T_cool_new - constants.T_cool_SS()))
            #FPI has converged, so the end time reactivity needs to be rho_0 for the next time step
            #and the end time power and DNP concentration needs to be the initial
            # conditions for the next time step
            rho_0 = rho_1
            r = np.append(r,rho_1)
            r_ext = np.append(r_ext,rho_ext)
            zeta = zeta_new
            P = P_new
            T_cool = T_cool_new
            T_fuel = T_fuel_new
            staggered_it = np.append(staggered_it,num_its)
            power_staggered = np.append(power_staggered,P_new)
        return power_staggered, r, r_ext, staggered_it