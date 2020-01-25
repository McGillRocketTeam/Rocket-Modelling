# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 16:20:40 2020

@author: sunge
"""
import CoolProp as cp
import math

from sympy import Symbol
from sympy.solvers import solve

class Nozzle:
    throat_area = 1
    throat_diameter = 1
    pressure = 1
    temperature = 1
    mass_flow_rate = 1

    def __init__(self, a, d, p, T):
        self.throat_area = a
        self.throat_diameter = d
        self.pressure = p
        self.temperature = T
    """
    4 methods below are based on Prof. Andrew Higgins' Fluids 2 notes
    """





    """
    get_outlet_mach: determines the outlet mach number based on the area expansion ratio
    #
    Takes inputs:
    outlet_area: outlet area
    throat_area: throat area
    gamma: heat capacity ratio (Cp/Cv), 1.4 for diatomic gas
    #
    Returns:
    outlet mach number
    """

    def get_outlet_mach(self, outlet_area, throat_area, gamma = 1.4):
        m = Symbol('M')
        sol = solve((outlet_area / throat_area) - (1 / m) * ((2 / (gamma + 1)) * (1 + (gamma - 1) / 2) * m ** 2) ** ((gamma + 1) / (2 * (gamma - 1))), m)
        print("Outlet mach number:",sol)
        return sol

    """
    get_mass_flow_rate: determines the mass flow rate based on combustion chamber pressure and temperature as well as
    throat area assuming sonic condition
    #
    Takes inputs:
    pressure: combustion chamber pressure
    temperature: combustion chamber temperature
    throat_area: throat area
    M_prime: throat mach number, will always be 1 as we assume sonic condition
    gamma: heat capacity ratio (Cp/Cv), 1.4 for diatomic gas
    R_prime = gas constant
    #
    Returns:
    mass_flow_rate: mass flow rate
    #
    """

    def get_mass_flow_rate(self, pressure, temperature, throat_area, M_prime = 1, gamma=1.4, R_prime = 8.314):
        m = Symbol('m')
        sol = solve(m - throat_area * (pressure/temperature) * (gamma / (R_prime / M_prime)) ** (0.5)
                    * (M_prime / (1 + ((gamma - 1) * M_prime ** 2) / 2) ** ((gamma + 1) / (2 * (gamma - 1)))))
        print("Mass flow rate:",sol)
        return sol

    """
    get_exit_pressure: determines the outlet pressure based on stagnation conditions (combustion chamber)
    #
    Takes inputs:
    p0: combustion chamber pressure (stagnation pressure)
    mach_num: outlet mach number
    gamma: heat capacity ratio (Cp/Cv), 1.4 for diatomic gas
    #
    Returns: 
    p = outlet pressure
    """

    def get_exit_pressure(self, p0,mach_num, gamma = 1.4):
        p = Symbol('p')
        sol = solve(p0/p - ((1+(gamma-1)* mach_num**2)/2) ** (gamma/(gamma-1)),m)
        print("Exit Pressure:",sol)
        return sol

    """
    get_thrust_from_nozzle: determines the thrust output from nozzle based on parameters obtained above
    #
    Takes inputs:
    m_dot =
    """

    def get_thrust_from_nozzle(self, m_dot, r_cnst, temp_not, p_e, p_o, p_amb, a_e, gamma=1.4):
        thrust_from_nozzle = m_dot * math.sqrt(((2 * gamma) / (gamma - 1)) * r_cnst * temp_not * (1 - (p_e / p_o) ** ((gamma - 1) / gamma))) + (p_e - p_amb) * a_e
        return thrust_from_nozzle

    """
    Determines the thrust output based on following input properties:
    #
    Takes inputs:
    throat_area: throat area
    outlet_area: outlet area
    pressure: combustion chamber pressure (Pstag)
    temperature: combustion chamber temperature (Tstag)
    #
    Return:
    thrust: thrust output of nozzle
    """
    def converge(self, throat_area, outlet_area, pressure, temperature,):
        mass_flow_rate = self.get_mass_flow_rate(pressure, temperature, throat_area)
        mach_num = self.get_outlet_mach(outlet_area,throat_area)
        exit_pressure = self.get_exit_pressure(pressure, mach_num)
        R = 8.314
        P_amb = 101324 #Pa
        thrust = self.get_thrust_from_nozzle(mass_flow_rate,R,temperature,exit_pressure,pressure,P_amb,outlet_area)
        print("Thrust (N):", thrust)
        return thrust