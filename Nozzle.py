# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 16:20:40 2020

@author: sunge
"""
import CoolProp as cp
import math

from sympy import Symbol
from sympy.solvers import solve
"""
Determines the nozzle outlet area or/and thrust output based on input variables
#
Inputs
gamma: heat capacity ratio (Cp/Cv) or (Cv + R / Cv)
throat_area: throat area
M_prime = throat mach number (assume sonic)
R_prime = universal gas constant 
mach_num = outlet mach number

"""


    def nozzle_area_ratio(self, outlet_area, throat_area, gamma):
        m = Symbol('M')
        sol = solve((outlet_area / throat_area) - (1 / m) * ((2 / (gamma + 1)) * (1 + (gamma - 1) / 2) * m ** 2) ** ((gamma + 1) / (2 * (gamma - 1))), m)
        print("Nozzle area ratio:",sol)
        return sol
    
    def mass_flow_rate(self, pressure, temperature, throat_area, M_prime = 1, gamma = 1.4, R_prime = 8.314):
        # find the mass flow rate of nozzle (assuming isentropic flow)
        # input = inlet pressure, inlet area
        # maximum mass flow rate obtained when inlet Mach num = 1
        m = Symbol('m')
        sol = solve(m - throat_area * (pressure/temperature) * (gamma / (R_prime / M_prime)) ** (0.5)
                    * (M_prime / (1 + ((gamma - 1) * M_prime ** 2) / 2) ** ((gamma + 1) / (2 * (gamma - 1)))))
        print("Mass flow rate:",sol)
        return sol
    
    def get_exit_pressure(p0,mach_num, gamma = 1.4):
        p = Symbol('p')
        sol = solve(p0/p - ((1+(gamma-1)* mach_num**2)/2) ** (gamma/(gamma-1)),m)
        print("Exit Pressure:",sol)
        return sol

    def get_thrust_from_nozzle(self, m_dot, r_cnst, temp_not, p_e, p_o, p_amb, a_e, gamma=1.4):
        thrust_from_nozzle = m_dot * math.sqrt(((2 * gamma) / (gamma - 1)) * r_cnst * temp_not * (1 - (p_e / p_o) ** ((gamma - 1) / gamma))) + (p_e - p_amb) * a_e
        return thrust_from_nozzle

