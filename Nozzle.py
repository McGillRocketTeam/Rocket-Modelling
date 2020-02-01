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
    A_throat = 1
    P_cc = 1
    T_cc = 1
    m_dot = 1

    def __init__(self, a, p, T):
        self.A_throat = a
        self.P_cc = p
        self.T_cc = T

    """
    4 methods below are based on Prof. Andrew Higgins' Fluids 2 notes
    """

    """
    get_outlet_mach: determines the outlet mach number based on the area expansion ratio
    #
    Takes inputs:
    A_exit: exit area [m^2]
    A_throat: throat area [m^2]
    gamma: heat capacity ratio (Cp/Cv), 1.4 for diatomic gas
    #
    Returns:
    outlet mach number
    """

    def get_outlet_mach(self, A_exit, A_throat, gamma = 1.4):
        m = Symbol('M')
        sol = solve((A_exit / A_throat) - (1 / m) * ((2 / (gamma + 1)) * ((1 + (gamma - 1) / 2) * m ** 2)) ** ((gamma + 1) / (2 * (gamma - 1))), m)
        M_out = float(sol[0])
        return M_out

    """
    get_mass_flow_rate: determines the mass flow rate based on combustion chamber P_cc and T_cc as well as
    throat area assuming sonic condition
    #
    Takes inputs:
    P_cc: combustion chamber pressure [Pa]
    T_cc: combustion chamber temperature [Pa]
    A_throat: throat area [m^2]
    M_prime: throat mach number, will always be 1 as we assume sonic condition
    gamma: heat capacity ratio (Cp/Cv), 1.4 for diatomic gas
    R = gas constant of N2O / Paraffin mixture at design O/F ratio [J/kgK]
    #
    Returns:
    m_dot: mass flow rate [kg/s]
    #
    """

    def get_mass_flow_rate(self, P_cc, T_cc, A_throat, R = 309.0878296654275, M_prime = 1, gamma=1.4):
        m = Symbol('m')
        sol = solve(m - A_throat * (P_cc / T_cc) * (gamma / (R / M_prime)) ** 0.5* (M_prime / (1 + ((gamma - 1) * M_prime ** 2) / 2) ** ((gamma + 1) / (2 * (gamma - 1)))))
        mdot = float(sol[0])
        return mdot
    """
    get_exit_pressure: determines the outlet P_cc based on stagnation conditions (combustion chamber)
    #
    Takes inputs:
    P_cc: combustion chamber pressure (stagnation pressure) [Pa]
    M: outlet mach number
    gamma: heat capacity ratio (Cp/Cv), 1.4 for diatomic gas
    #
    Returns: 
    P_exit: exit pressure [Pa]
    """

    def get_exit_pressure(self, P_cc, M, gamma = 1.4):
        p = Symbol('p')
        sol = solve(P_cc / p - (1 + ((gamma - 1) * M ** 2) / 2) ** (gamma / (gamma - 1)), p)
        P_exit = float(sol[0])
        return P_exit

    """
    get_thrust_from_nozzle: determines the thrust output from nozzle based on parameters obtained above
    #
    Takes inputs:
    m_dot: mass flow rate [kg/s]
    T_cc: combustion chamber temperature [K]
    P_exit: exit pressure [Pa]
    P_cc: combustion chamber temperature [Pa]
    P_amb: ambient pressure [Pa]
    A_exit: exit area [m^2]
    R = gas constant of N2O / Paraffin mixture at design O/F ratio [J/kgK]
    gamma: heat capacity ratio (Cp/Cv), 1.4 for diatomic gas
    #
    Returns:
    thrust: thrust output [N]
    """

    def get_thrust_from_nozzle(self, m_dot, T_cc, P_exit, P_cc, P_amb, A_exit, R = 309.0878296654275, gamma=1.4):
        thrust_from_nozzle = m_dot *(((2 * gamma) / (gamma - 1)) * R * T_cc * (1 - (P_exit / P_cc) ** ((gamma - 1) / gamma)))**0.5 + (P_exit - P_amb) * A_exit
        return thrust_from_nozzle

    """
    Determines the thrust output based on following input properties:
    #
    Takes inputs:
    A_throat: throat area [m^2]
    A_exit: outlet area [m^2]
    P_cc: combustion chamber pressure (Pstag) [Pa]
    T_cc: combustion chamber temperature (Tstag) [Pa]
    P_amb: ambient pressure [Pa]
    #
    Return:
    thrust: thrust output of nozzle [N]
    """

    def converge(self, A_throat, A_exit, P_cc, T_cc, P_amb):
        m_dot = self.get_mass_flow_rate(P_cc, T_cc, A_throat)
        print(m_dot)
        mach_num = self.get_outlet_mach(A_exit, A_throat)
        print(mach_num)
        exit_pressure = self.get_exit_pressure(P_cc, mach_num)
        print(exit_pressure)
        thrust = self.get_thrust_from_nozzle(m_dot, T_cc, exit_pressure, P_cc, P_amb, A_exit)
        print("Thrust [N]:", thrust)
        return thrust


"""
TEST SCRIPT FROM BELOW

nozzleTest = Nozzle(2,300000,300)
testConverge = nozzleTest.converge(nozzleTest.A_throat,4,nozzleTest.P_cc, nozzleTest.T_cc,101325)

OUTPUT
77.894953578635
1.148698354997035
132167.55424695372
Thrust [N]: 152046.59462488594
"""