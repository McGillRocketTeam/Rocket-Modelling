# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 16:13:36 2020

@author: sunge
"""
from math import sqrt
import CoolProp as cp

class Injector:
    area = 1
    length = 1
    mass = 1
    pressure = 1
    mass_flow_rate = 1
    temp = 1

    def __init__ (self,A,L,M):
        self.area = A_o
        self.length = L
        self.mass = M

    def design (self, P_tank, rho_tank, h_tank,  P_cc, rho_cc, h_cc, P_vap, C_d, m_dot):
        #break up parts of the equations for readability
        d_P = P_tank-P_cc
        d_h = h_tank-h_cc

        kappa = sqrt((d_P) / (P_vap - P_cc))
        
        m_inc_coef = (kappa/(1+kappa))*sqrt(2*rho_tank*d_P)
        m_hem_coef = (1/(1+kappa))*sqrt(2*rho_cc*d_h)

        area = sqrt(m_dot/((C_d**2)*(m_inc_coef + m_hem_coef))
        return area


    def converge (self, P_tank, rho_tank, h_tank,  P_cc, rho_cc, h_cc,  P_vap, C_d):
        #break up parts of the equations for readability
        d_P = P_tank - P_cc
        d_h = h_tank - h_cc

        kappa = sqrt(d_P/ (P_vap - P_cc))

        m_inc_coef = (kappa / (1 + kappa)) * C_d* self.area * sqrt(2 * rho_tank * d_P)
        m_hem_coef = (1/(1+kappa))*C_d*self.area*sqrt(2*rho_cc*d_h)

        m_dot = C_d*self.area*sqrt(m_inc_coef+m_hem_coef)
        return m_dot
        
