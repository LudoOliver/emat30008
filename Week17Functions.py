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
    u = np.array(u,dtype=np.float64)
    du1 = Beta*u[0]-u[1]-u[0]*(u[0]**2+u[1]**2)
    du2 = u[0]+Beta*u[1]-u[1]*(u[0]**2+u[1]**2)
    return np.array([du1,du2],dtype=np.float64)

def ModifiedBetaFormHopf(u,t,Beta=1):
    #u = np.array(u,dtype=np.float64)

    #try:
        #du1 = Beta*u[0]-u[1]+u[0]*(u[0]**2+u[1]**2)*(1-(u[0]**2+u[1]**2))
        #du1
    #except:
        #print(Beta*u[0]-u[1]+u[]*(u[0]**2+u[1]**2)*(1-(u[0]**2+u[1]**2)))
    #print(type(du1))
    #if type(du1)!=np.float64:
    #   print("error")
    #    return
    #exponent1=(((u[0])**2)+(u[1])**2)
    #if np.isnan(exponent1):
    #    print(u[0])
    #   print(u[0**2])
    #    print(u[1**2])
    #exponent2= exponent1**2
    #if np.isnan(exponent2):
        #print(exponent1)
        #return
    du1 = Beta*u[0]-u[1]+u[0]*((u[0])**2+(u[1])**2)-u[0]*(((u[0])**2+(u[1])**2)**2)
    du2 = u[0]-Beta*u[1]+u[1]*((u[0])**2+(u[1])**2)-u[1]*(((u[0])**2+(u[1])**2)**2)
    return np.array([du1,du2])#,dtype=np.float64)

