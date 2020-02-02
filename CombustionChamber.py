# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 16:15:22 2020

This piece of code creates a CombustionChamber class, which tracks fuel mass
and geometry over time, as well as having a method to determine the regression
based on experimental constants.

The model here is based on the 1D model of regression, which assumes:
    1. no space dependency for port regression
    2. No effect of pressure on regression rate

See Cantwell's chapter on hybrid rocket engines for more info

@author: jotisl, sunge
"""
import NumPy as np
class CombustionChamber:
    outer_radius = 1
    inner_radius = 1
    length = 1
    pressure = 1
    temp = 1
    mass_flow_rate = 1

    # two constants: defined in diff. file or from literature
    a_ballistic = 0.0000927
    n_ballistic = 0.62

    def __init__(self, a, b, c, d, e, f):
        self.outer_radius = a
        self.inner_radius = b
        self.grain_length = c
        self.pressure = d
        self.temperature = e
        self.rho_fuel = f
        self.m_fuel = self.rho_fuel*self.grain_length*np.pi*(self.outer_radius**2 - self.inner_radius**2)

    def chamberTemp(self, OF_ratio, P_cc):
        #here, call CEA to extract combustion equilibrium stuff
        return self.temperature

    def update(self, dt, m_dot_fuel):
        m_lost = m_dot_fuel*dt
        self.m_fuel -= m_lost
        delta_vol = m_lost*self.rho_fuel
        delta_r = delta_vol / (4*np.pi*self.inner_radius**2)
        self.inner_radius += delta_r #here the regression number is positive in the r-direction

    def converge(self, m_dot_ox, P_cc):
        r_dot_fuel = self.n_ballistic*(m_dot_ox/(np.pi*self.inner_radius**2))**self.n_ballistic
        m_dot_fuel = self.rho_fuel*self.grain_length*np.pi*2*self.inner_radius*r_dot_fuel
        
        OF_ratio = m_dot_ox/m_dot_fuel
        OF_ratio_f = (m_dot_ox**(1-self.n_ballistic)*(2*self.inner_radius)**(2*self.n_ballistic-1))/(4**self.n_ballistic*np.pi**(1-self.n_ballistic)*self.n_ballistic*self.rho_fuel*self.grain_length)
        
        self.temperature = self.chamberTemp(OF_ratio, P_cc)
        
        
        return m_dot_fuel, self.temperature
        
        