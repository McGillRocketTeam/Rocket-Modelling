# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 16:06:20 2020

@author: sunge, crw
"""

from CoolProp.CoolProp import PropsSI
import numpy as np

class Tank:

    DEBUG_VERBOSITY = 1

    total_volume = 1
    oxidizer_mass = 1
    initial_temp = 1
    tank_dry_mass = 1
    outlet_diameter = 1

    def __init__(self):
        self.V_tot = 0.0232
        self.m_ox = 10
        self.v_specific = self.V_tot/self.m_ox

        self.tank_id = .1346
        self.T_tank = 330  # Initial temperature, kelivin
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
        T = PropsSI("T", "H", self.h, "D", 1 / self.v_specific, 'N2O')

        if self.DEBUG_VERBOSITY > 0:
            print("***DEBUG*** [Tank.update] Tank Oxidiser Mass = ", self.m_ox)

        if self.DEBUG_VERBOSITY > 1:
            print("***DEBUG*** [Tank.update] Tank Temperature = ", T)
            p = PropsSI("P", "T", self.T_tank, "D", 1 / self.v_specific, 'N2O')
            print("***DEBUG*** [Tank.update] Tank Pressure = ", p)

        return T, rho
