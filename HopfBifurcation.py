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
    #try: 
        #print(du1,du2,'*')
    #except:
        #print("Error")
        #return 
    return np.array([du1,du2])
# def ShootingHopfTest():
#     FoundX,FoundPeriod= Week16General.Shooting(NormalFormHopf, X0, T0, StepSize=0.0001)
#     Error = math.sqrt(Beta)*math.sin(math.acos(FoundX[0])/math.sqrt(Beta))-FoundX[1]
#     if abs(Error)<0.05:
#         print("Success")
#         return True
#     else:
#         print("Failiure")
#         return False 
X0 =[0.1,0.1]
T0 = 8
Tolerance =0.05
Xn,Tn= Week16General.Shooting(NormalFormHopf, X0, T0, StepSize=0.001)
SolnXArray,SolnTArray = ODESolver.Solve_to(NormalFormHopf,Xn,[0,Tn],0.001,ODESolver.RungeKutta4)
TrueXArray = np.vstack((math.sqrt(Beta)*np.cos(SolnTArray-math.pi),math.sqrt(Beta)*np.sin(SolnTArray-math.pi)))
ComparisonArray = np.isclose(SolnXArray,TrueXArray)#,rtol=1e-3,atol=1e-1)

            

#%%
def ShootingHopfTestProper(Tolerance=0.05):
    X0 =[0.1,0.1]
    T0 = 8
    Xn,Tn= Week16General.Shooting(NormalFormHopf, X0, T0, StepSize=0.001)
    SolnXArray,SolnTArray = ODESolver.Solve_to(NormalFormHopf,Xn,Tn,0.001,ODESolver.RungeKutta4)
    TrueXArray = np.vstack((math.sqrt(Beta)*np.cos(SolnTArray+math.pi/2),math.sqrt(Beta)*np.sin(SolnTArray+math.pi/2)))
    ComparisonArray = np.isclose(SolnXArray,TrueXArray,atol=Tolerance)
    print(ComparisonArray)
    if np.all(ComparisonArray):
        print("Success")
        return 
    else:
        print("Failiure")
        return    
        
Answer=ShootingHopfTestProper(0.05)        
#%%   
X0 =[0.1,0.1]
T0 = 8

FoundX,FoundPeriod= Week16General.Shooting(NormalFormHopf, X0, T0, StepSize=0.0001)
#%%
#= X,T
#print(Solution)
#FoundICs = [FoundX,FoundY]
FoundTime = [0, FoundPeriod]
x,t = ODESolver.Solve_to(NormalFormHopf,FoundX,FoundTime,0.00001,ODESolver.RungeKutta4)
plt.plot(x[0,:],x[1,:]) 
plt.figure()
#plt.plot(t,x[1,:]) #Shows that period is roughly 32 for b= 0.1 and a 
#plt.show()