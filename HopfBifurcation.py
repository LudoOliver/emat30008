# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 13:35:16 2024

@author: xd21736
"""

import scipy
import math
import numpy as np
import matplotlib.pyplot as plt
import ODESolver
import Week16General

Beta=0.9
Sigma=-1


def NormalFormHopf(u,t):
    du1 = Beta*u[0]-u[1]+Sigma*u[0]*(u[0]**2+u[1]**2)
    du2 = u[0]+Beta*u[1]+Sigma*u[1]*(u[0]**2+u[1]**2)
    return np.array([du1,du2])

def ThreeDimFormHopf(u,t):
    du1 = Beta*u[0]-u[1]+Sigma*u[0]*(u[0]**2+u[1]**2)
    du2 = u[0]+Beta*u[1]+Sigma*u[1]*(u[0]**2+u[1]**2)
    du3 = -u[2]
    return np.array([du1,du2,du3])

X0 =[0.1,0.1]
T0 = 8
Tolerance =0.05


#%%
def ShootingHopfTestProper(Tolerance=0.05,X0 =[0.1,0.1],T0=8):
    """

    Parameters
    ----------
    Tolerance : Scalar, optional
        DESCRIPTION. The default is 0.05.
    X0 : , optional
        DESCRIPTION. The default is [0.1,0.1].
    T0 : TYPE, optional
        DESCRIPTION. The default is 8.

    Returns
    -------
    State : Boolean value of test outcome

    """
    Xn,Tn= Week16General.Shooting(NormalFormHopf, X0, T0, StepSize=0.001)
    SolnXArray,SolnTArray = ODESolver.Solve_to(NormalFormHopf,Xn,[0,Tn],0.001,ODESolver.RungeKutta4)
    TrueXArray = np.vstack((math.sqrt(Beta)*np.cos(SolnTArray-math.pi),math.sqrt(Beta)*np.sin(SolnTArray-math.pi)))
    #-pi as phase condtion requires du0 to be 0
    ComparisonArray = np.isclose(SolnXArray,TrueXArray,atol=Tolerance)
    if np.all(ComparisonArray):
        print("Solution is within tolerance of the real limit cycle")
        return 1
    else:
        print("Out of Tolerance")
        return 0  
    
def ShootingThreeHopf(Tolerance=0.5,X0=[0.1,0.1,1e-9],T0=8):
    #print(Tolerance,X0,T0)
    Xn,Tn= Week16General.Shooting(ThreeDimFormHopf, X0, T0, StepSize=0.001)
    print(Xn)
    SolnXArray,SolnTArray = ODESolver.Solve_to(ThreeDimFormHopf,Xn,[0,Tn],0.001,ODESolver.RungeKutta4)
    TrueXArray = np.vstack((math.sqrt(Beta)*np.cos(SolnTArray-math.pi),math.sqrt(Beta)*np.sin(SolnTArray-math.pi),Xn[-1]*np.exp(-SolnTArray)))
    #-pi as phase condtion requires du0 to be 0
    ComparisonArray = np.isclose(SolnXArray,TrueXArray,atol=Tolerance)

    if np.all(ComparisonArray):
        print("Solution is within tolerance of the real limit cycle")
        return 1
    else:
        print("Out of Tolerance")
        return 0  
    
        
Answer=ShootingHopfTestProper(0.05)      
Answer2=ShootingThreeHopf(0.1)  
#%%   
X0 =[0.1,0.1]
T0 = 8
#%%



#plt.plot(t,x[1,:]) #Shows that period is roughly 32 for b= 0.1 and a 
#plt.show()