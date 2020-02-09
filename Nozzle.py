# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 16:20:40 2020

@author: sunge, jyoo, crw
"""
import CoolProp as cp
import math

from sympy import Symbol
from sympy.solvers import solve


class Nozzle:

    DEBUG_VERBOSITY = 0


    A_throat = 1
    P_cc = 1
    T_post_comb = 1
    m_dot = 1

    def __init__(self):
        self.A_throat = 0.5
        self.A_exit = 0.1
        #self.P_cc = p
        #self.T_post_comb = T
        self.gamma = 1.4
        self.gas_constant = 8.314*1e-3/27.4
        self.thrust=0

    """
    4 methods below are based on Prof. Andrew Higgins' Fluids 2 notes
    """

    """
    get_outlet_mach: determines the outlet mach number based on the area expansion ratio
    #
    Takes inputs:
    gamma: heat capacity ratio (Cp/Cv), 1.4 for diatomic gas
    #
    Returns:
    outlet mach number
    """

    def get_outlet_mach(self):
        m = Symbol('M')
        A_throat = self.A_throat
        A_exit = self.A_exit
        gamma = self.gamma
        sol = solve((A_exit / A_throat) - (1 / m) * ((2 / (gamma + 1)) * ((1 + (gamma - 1) / 2) * m ** 2)) ** ((gamma + 1) / (2 * (gamma - 1))), m)
        M_out = float(sol[0])
        return M_out

    """
    get_mass_flow_rate: determines the mass flow rate based on combustion chamber P_cc and T_post_comb as well as
    throat area assuming sonic condition
    #
    Takes inputs:
    P_cc: combustion chamber pressure [Pa]
    T_post_comb: combustion chamber temperature [Pa]
    R = gas constant of N2O / Paraffin mixture at design O/F ratio [J/kgK]
    #
    Returns:
    m_dot: mass flow rate [kg/s]
    #
    """

    def get_mass_flow_rate(self, P_cc, T_post_comb):
        m = Symbol('m')
        M_prime = 1
        A_throat = self.A_throat
        R = self.gas_constant
        gamma = self.gamma
        sol = solve(m - A_throat * (P_cc / T_post_comb) * (gamma / (R / M_prime)) ** 0.5* (M_prime / (1 + ((gamma - 1) * M_prime ** 2) / 2) ** ((gamma + 1) / (2 * (gamma - 1)))))
        mdot = float(sol[0])

        if self.DEBUG_VERBOSITY > 0:
            print("[Nozzle.get_mass_flow_rate] mdot = ", mdot)

        return mdot
    """
    get_exit_pressure: determines the outlet P_cc based on stagnation conditions (combustion chamber)
    #
    Takes inputs:
    P_cc: combustion chamber pressure (stagnation pressure) [Pa]
    M: outlet mach number
    #
    Returns: 
    P_exit: exit pressure [Pa]
    """

    def get_exit_pressure(self, P_cc, M):
        p = Symbol('p')
        gamma = self.gamma
        sol = solve(P_cc / p - (1 + ((gamma - 1) * M ** 2) / 2) ** (gamma / (gamma - 1)), p)
        P_exit = float(sol[0])

        if self.DEBUG_VERBOSITY > 1:
            print("[Nozzle.get_exit_pressure] p_exit = ", P_exit)

        return P_exit

    """
    get_thrust_from_nozzle: determines the thrust output from nozzle based on parameters obtained above
    #
    Takes inputs:
    m_dot: mass flow rate [kg/s]
    T_post_comb: combustion chamber temperature [K]
    P_exit: exit pressure [Pa]
    P_cc: combustion chamber temperature [Pa]
    P_amb: ambient pressure [Pa]
    #
    Returns:
    thrust: thrust output [N]
    """

    def get_thrust_from_nozzle(self, m_dot, T_post_comb, P_exit, P_cc, P_amb):
        R = self.gas_constant
        gamma = self.gamma
        A_exit = self.A_exit
        thrust_from_nozzle = m_dot *(((2 * gamma) / (gamma - 1)) * R * T_post_comb * (1 - (P_exit / P_cc) ** ((gamma - 1) / gamma)))**0.5 + (P_exit - P_amb) * A_exit
        self.thrust = thrust_from_nozzle

        if self.DEBUG_VERBOSITY > 0:
            print("***DEBUG*** [Nozzle.get_thrust_from_nozzle] Thrust = ", thrust_from_nozzle)

        return thrust_from_nozzle

    """
    Determines the thrust output based on following input properties:
    #
    Takes inputs:
    A_throat: throat area [m^2]
    A_exit: outlet area [m^2]
    P_cc: combustion chamber pressure (Pstag) [Pa]
    T_post_comb: combustion chamber temperature (Tstag) [Pa]
    P_amb: ambient pressure [Pa]
    #
    Return:
    thrust: thrust output of nozzle [N]
    """

    def converge(self, m_dot_actual, T_post_comb, P_cc):
        m_dot_choke = self.get_mass_flow_rate(P_cc, T_post_comb)
        return m_dot_choke

    def update(self, dt):
        pass
"""
TEST SCRIPT FROM BELOW

nozzleTest = Nozzle(2,300000,300)
testConverge = nozzleTest.converge(nozzleTest.A_throat,4,nozzleTest.P_cc, nozzleTest.T_post_comb,101325)

OUTPUT
77.894953578635
1.148698354997035
132167.55424695372
Thrust [N]: 152046.59462488594
"""
