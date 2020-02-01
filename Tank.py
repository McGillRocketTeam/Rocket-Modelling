# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 16:06:20 2020

@author: sunge
"""

from CoolProp.CoolProp import PropsSI
import Numpy as np

class Tank:
    total_volume = 1
    oxidizer_mass = 1
    initial_temp = 1
    tank_dry_mass = 1
    outlet_diameter = 1

    def __init__(self, init):
        self.V_tot = a
        self.m_ox = b
        self.v_specific = self.V_tot/self.m_ox

        self.tank_id = c
        self.T_tank = d  # Initial temperature
        self.h = PropsSI("H", "T", self.T_tank, "D", 1/self.v_specific, 'N2O')

# example usage of PropsSI: PropsSI("T", "P", 101325, "Q", 0, "Water")
# finds temperature of water at pressure 101325 Pa and Quality 0 (which is equivalent to boiling temp)

    def converge(self):
        pass

    def update(self, dt, m_dot_ox):
        # used to move timestep forward
        self.m_ox -= m_dot_ox * dt
        self.v_specific = self.V_tot / self.m_ox
        rho = 1 / self.v_specific
        T = PropsSI("T", "H", self.T_tank, "D", 1 / self.v_specific, 'N2O')
        return T, rho
