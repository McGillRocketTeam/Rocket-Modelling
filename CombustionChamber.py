# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 16:15:22 2020

@author: sunge
"""
import NumPy as np
class CombustionChamber:
    outer_radius = 1
    inner_radius = 1
    length = 1
    pressure = 1
    temp = 1
    mass_flow_rate = 1

    # two constant: defined in diff. file or
    a_ballistic = 1
    n_ballistic = 1

    def __init__(self, a, b, c, d):
        self.outer_radius = a
        self.inner_radius = b
        self.grain_length = q
        self.pressure = c
        self.temperature = d
        self.rho_fuel = f
        self.m_fuel = self.rho_fuel*self.grain_length*np.pi*(self.outer_radius**2 - self.inner_radius**2)

    def update(self, dt, m_dot_fuel):
        m_lost = m_dot_fuel*dt
        self.m_fuel -= m_lost
        delta_vol = m_lost*self.rho_fuel
        delta_r = delta_vol / (4*np.pi*self.inner_radius**2)
        self.inner_radius -= delta_r

