# -*- coding: utf-8 -*-
"""
Created on Sat Jan 18 13:26:42 2020

@author: sunge
"""
import CoolProp.CoolProp as CP
fluid = 'Water'
pressure_at_critical_point = CP.PropsSI(fluid,'pcrit')
# Massic volume (in m^3/kg) is the inverse of density
# (or volumic mass in kg/m^3). Let's compute the massic volume of liquid
# at 1bar (1e5 Pa) of pressure
vL = 1/CP.PropsSI('D','P',1e5,'Q',0,fluid)
# Same for saturated vapor
vG = 1/CP.PropsSI('D','P',1e5,'Q',1,fluid)