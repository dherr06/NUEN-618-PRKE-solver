class ragusa_code:
    ################################################################################
    # Fuel thermal conductivity function
    
    def k_fuel(temp):
        y = 1.5443e-6 * temp**2 - 0.0050449 * temp + 6.3396 # [W/m C]
        dy = 1.5443e-6 * 2 * temp - 0.0050449;
        return y, dy

    ################################################################################
    # Fuel density times heat capacity function
    
    def rhocp_fuel(temp):
        rho = 10412.0 # [kg/m^3]
        cp = 340.0 # [J/kg C]
        y = rho * cp
        dy = 0
        return y, dy

    ################################################################################
    # Moderator heat capacity function

    def cp_mod(x):
        y = 6.736093e-7 * x**4 - 8.224535e-4 * x**3 + 3.767976e-1 * x**2 - \
            7.673778e1 * x + 5.865020e3
        dy = 6.736093E-07 * 4 * x**3 - 8.224535E-04 * 3 * x**2 + \
             3.767976E-01 * 2 * x - 7.673778E+01
        y *= 1e3
        dy *= 1e3
        return y, dy

    ################################################################################
    # Moderator density function
    
    def rho_mod(temp):
        y = -0.0196 * temp**2 + 9.7317 * temp - 430.65 # [kg/m^3]
        dy = -0.0196 * 2 * temp + 9.7317
        return y, dy
    
    ################################################################################
    # Moderator density times heat capacity
    
    def rhocp_mod(temp):

        def rho_mod(temp):
            y = -0.0196 * temp**2 + 9.7317 * temp - 430.65 # [kg/m^3]
            dy = -0.0196 * 2 * temp + 9.7317
            return y, dy

        def cp_mod(x):
            y = 6.736093e-7 * x**4 - 8.224535e-4 * x**3 + 3.767976e-1 * x**2 - \
                7.673778e1 * x + 5.865020e3
            dy = 6.736093E-07 * 4 * x**3 - 8.224535E-04 * 3 * x**2 + \
                 3.767976E-01 * 2 * x - 7.673778E+01
            y *= 1e3
            dy *= 1e3
            return y, dy

        rho, drho = rho_mod(temp)
        cp, dcp = cp_mod(temp)
        y = rho * cp
        dy = drho * cp + rho * dcp
        return y, dy
    