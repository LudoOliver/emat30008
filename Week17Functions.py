# -*- coding: utf-8 -*-
"""
Week17 Functions
Created on Mon Feb 26 14:26:31 2024

@author: farka
"""
from CommonModules import *

def Cubic(x,t,c=1):
    return x**3-x+c

def BetaFormHopf(u,t,Beta=1):
    du1 = Beta*u[0]-u[1]-u[0]*(u[0]**2+u[1]**2)
    du2 = u[0]+Beta*u[1]-u[1]*(u[0]**2+u[1]**2)
    return np.array([du1,du2])

def ModifiedBetaFormHopf(u,t,Beta=1):
    du1 = Beta*u[0]-u[1]+u[0]*(u[0]**2+u[1]**2)*(1-(u[0]**2+u[1]**2))
    du2 = u[0]+Beta*u[1]-u[1]*(u[0]**2+u[1]**2)*(1-(u[0]**2+u[1]**2))
    return np.array([du1,du2])

