# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 16:06:20 2020

@author: sunge, crw
"""

from CoolProp.CoolProp import PropsSI, get_phase_index, PhaseSI
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
        #self.tank_id = .1346
        self.T_tank = 35 + 273  # Initial temperature, kelvin

        self.v_specific_mix = self.V_tot/self.m_ox
        self.H = self.m_ox*PropsSI("H", "T", self.T_tank, "D", 1/self.v_specific_mix, 'N2O')
        self.rho_liquid = PropsSI("D", "T", self.T_tank, "Q", 0, "N2O")

# example usage of PropsSI: PropsSI("T", "P", 101325, "Q", 0, "Water")
# finds temperature of water at pressure 101325 Pa and Quality 0 (which is equivalent to boiling temp)

    def converge(self):
        pass

    def update(self, dt, m_dot_ox):
        # used to move timestep forward

        # MODELLING ASSUMPTIONS:
        # only liquid phase is removed from tank, so long as we aren't critical
        # the mixture enthalpy changes as func of time bc of this (h_vap > h_liq)
        # Assume no other enthalpy sinks though (e.g. adiabatic)
        # Simulation ends when no longer two-phase (e.g. gas, critical)

        # Check to ensure phase is twophase, otherwise end simulation
        current_phase_index = PropsSI("Phase", "T", self.T_tank, "D", 1/self.v_specific_mix, "N2O")
        twophase_phase_index = get_phase_index("phase_twophase")
        if current_phase_index == twophase_phase_index:

            # calculate specific enthalpy of liquid phase
            h_liq = PropsSI("H", "T", self.T_tank, "Q", 0, 'N2O')

            # Update mass and total enthalpy accordingly
            self.m_ox -= m_dot_ox * dt
            self.H -= (m_dot_ox * dt) * h_liq

            # Re-evaluate other properties
            self.v_specific_mix = self.V_tot / self.m_ox
            h_mix = self.H/self.m_ox
            self.T_tank = PropsSI("T", "H", h_mix, "D", 1 / self.v_specific_mix, 'N2O')
            self.rho_liquid = PropsSI("D", "T", self.T_tank, "Q", 0, "N2O")

            if self.DEBUG_VERBOSITY > 0:
                print("***DEBUG*** [Tank.update] Tank Oxidiser Mass = ", self.m_ox)

            if self.DEBUG_VERBOSITY > 1:
                print("***DEBUG*** [Tank.update] Tank Temperature = ", self.T_tank)
                print("***DEBUG*** [Tank.update] Tank Total Enthalpy = ", self.H)
                p = PropsSI("P", "T", self.T_tank, "D", 1 / self.v_specific_mix, 'N2O')
                print("***DEBUG*** [Tank.update] Tank Pressure = ", p)

        else:
            print("[Tank] Nitrous is no longer two-phase mixture, ending simulation")
            self.rho_liquid = -1
            # print(PhaseSI("T", self.T_tank, "D", 1 / self.v_specific_mix, "N2O"))

        return self.T_tank, self.rho_liquid
