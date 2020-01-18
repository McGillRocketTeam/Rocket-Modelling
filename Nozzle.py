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
        print(sol)