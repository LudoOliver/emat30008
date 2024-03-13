# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 20:26:26 2024

@author: Admin
"""

from CommonModules import *
import ODESolver

def FirstOrderFiniteDiff(Func,StepSize,LocationPoint):
    dx = (Func(LocationPoint+StepSize)-Func(LocationPoint-StepSize))
    dt = (2*StepSize)
    return dx/dt
def SecondOrderFiniteDiff(Func,StepSize,LocationPoint):
    ddx = Func(LocationPoint+StepSize)-2*Func(LocationPoint)+Func(LocationPoint-StepSize)
    dtt = StepSize**2
    return ddx/dtt

Alpha = 1
Beta = 5
Bounds = [Alpha,Beta]

def FiniteSolvePoisson(InputApprox,Bounds=Bounds,MaxStepSize=0.01):
    Gamma1 = 1 #definitely should change
    Gamma2 = 10
    NPoints = math.ceil((max(Bounds)-min(Bounds))/MaxStepSize)
    TrueStepSize = max(Bounds)-min(Bounds)/NPoints
    CoeffMatrix = np.zeros([NPoints,NPoints])
    #for i in range(LenCoeffMatrix):
    CoeffVector = np.array([-1,0,1])
    CoeffMatrix[0,(0,1)] = -2,1

    for i in range(1,len(CoeffMatrix)-1):
        CoeffMatrix[i,i+CoeffVector] = (1,-2,1)
        #print(i)
    CoeffMatrix[-1,(-2,-1)] = (1,-2)
    print(CoeffMatrix)
    
    U = np.zeros(NPoints)
    U[0] = Gamma1
    U[-1] = Gamma2
    
    def EqnToSolve(SolnArray):
        return CoeffMatrix.dot(SolnArray)+U
        
    Ans = scipy.optimize.root(EqnToSolve,U).x
    #Ans = Ans
    plt.figure()
    plt.plot(np.linspace(Gamma1,Gamma2,num=NPoints),Ans,label="mysoln")
    plt.plot(np.linspace(Gamma1,Gamma2),(Gamma2-Gamma1)*np.linspace(0,1)+Gamma1+0.1)
    plt.legend()
    plt.show()
    print(Ans)
    return 0

if __name__ == "__main__":
    FiniteSolvePoisson(1)