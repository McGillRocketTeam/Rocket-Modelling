# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 16:20:40 2020

@author: sunge
"""
from sympy import Symbol
from sympy.solvers import solve

class Nozzle:
    throat_area = 1
    outlet_area = 2
    mass = 3
    def __init__(self,a,b,c):
        self.throat_area = a
        self.outlet_area = b
        self.mass = c

    def flow_area_realtion(self, a, a_star, gamma):
        m = Symbol('M')
        sol = solve((a / a_star) - (1 / m) * ((2 / (gamma + 1)) * (1 + (gamma - 1) / 2) * m ** 2) ** ((gamma + 1) / (2 * (gamma - 1))), m)
        print("Flow area realton:",sol)
        return sol
    
    def mass_flow_rate(self, pressure, temperature,throat_area, mech_num = 1,
                       gamma = 1.4, R_prime = 8.314, M_prime = 30):
        m = Symbol('m')
        sol = solve(m - throat_area * (pressure/temperature) * (gamma/(R_prime/M_prime))**(0.5)
        * (mech_num/(1+((gamma-1)*mech_num**2)/2) ** ((gamma+1)/(2*(gamma-1)))))))))))
        print("Mass flow rate:",sol)
        return sol
        
        