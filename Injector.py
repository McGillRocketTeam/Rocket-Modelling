# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 16:13:36 2020

@author: sunge
"""
from math import sqrt
from CoolProp.CoolProp import PropsSI


class Injector:
    area = 1
    length = 1
    mass = 1
    pressure = 1
    mass_flow_rate = 1
    temp = 1

    def __init__ (self):
        self.area = 0.000048
        # self.length = L
        # self.mass = M
        self.C_d = 0.62

    """
    Determines the injection area based on input variables
    #
    Takes Input:
    T_tank: temperature of tank
    rho_tank: density of nitrous liquid vapour mixture
    T_cc: temperature of the combustion chamber
    rho_cc: density of propellent mixture in combustion chamber
    P_vap: vapour pressure of nitrous
    C_d: coefficient of
    m_dot: mass flow rate
    #
    Return:
    area: injection area: area of the little holes which inject fuel
    """
    def design(self, T_tank, rho_tank,  T_cc, rho_cc, P_vap, C_d, m_dot):
        # Find missing properties from CoolProp
        P_tank = PropsSI('P', 'T', T_tank, 'D', rho_tank, 'N2O')
        h_tank = PropsSI('H', 'T', T_tank, 'D', rho_tank, 'N2O')

        P_cc = PropsSI('P', 'T', T_cc, 'D', rho_cc, 'N2O')
        h_cc = PropsSI('H', 'T', T_cc, 'D', rho_cc, 'N2O')

        # break up parts of the equations for readability
        d_P = P_tank-P_cc
        d_h = h_tank-h_cc

        kappa = sqrt(d_P / (P_vap - P_cc))
        
        m_inc_coef = (kappa/(1+kappa))*sqrt(2*rho_tank*d_P)
        m_hem_coef = (1/(1+kappa))*sqrt(2*rho_cc*d_h)

        area = sqrt(m_dot/((C_d**2)*(m_inc_coef + m_hem_coef)))
        
        return area

    """
    determine the mass flow rate through the injector based on input properties
    #
    Takes Input:
    T_tank: temperature of tank
    rho_tank: density of nitrous liquid vapour mixture
    T_cc: temperature of the combustion chamber
    rho_cc: density of propellent mixture in combustion chamber
    P_vap: vapour pressure of nitrous
    C_d: coefficient of
    #
    Return:
    m_dot: mass flow rate through the injector
    """
    def converge(self, T_tank, rho_tank,  T_cc, rho_cc):
        # Find missing properties from CoolProp
        P_tank = PropsSI('P', 'T', T_tank, 'D', rho_tank, 'N2O')
        h_tank = PropsSI('H', 'T', T_tank, 'D', rho_tank, 'N2O')

        P_cc = PropsSI('P', 'T', T_cc, 'D', rho_cc, 'N2O')
        h_cc = PropsSI('H', 'T', T_cc, 'D', rho_cc, 'N2O')

        # break up parts of the equations for readability
        d_P = P_tank - P_cc
        d_h = h_tank - h_cc

        kappa = sqrt(d_P / (P_tank - P_cc))

        m_inc_coef = (kappa / (1 + kappa)) * self.C_d * self.area * sqrt(2 * rho_tank * d_P)
        m_hem_coef = (1/(1+kappa))*self.C_d*self.area*sqrt(2*rho_cc*d_h)

        m_dot = self.C_d*self.area*sqrt(m_inc_coef + m_hem_coef)
        return m_dot

    def update(self, dt):
        pass
