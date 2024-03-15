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

def NoSourceTerm(x,U):
    return 0

def SimpleSourceTerm(x,U):
    return 1 # Probably could just be 

def BratuTerm(x,U,mu=0.1):
    return np.exp(mu*U)

def FiniteSolvePoisson(Bounds,
                        DiffusionConstant=1,Reaction=NoSourceTerm,
                        NPoints=100,Guess=0,
                        Robin=0,
                        Neuman=0):
    
    
    # Just realised my step size idea sucks should just use npoints lol
    # Thats fine tho should be easy fix
    # Unsure about how different lengths should be
    # Should also be able to do different length scales some point
    
    X0 = Bounds[0,0]
    U0 = Bounds[1,0]
    XN = Bounds[0,1]
    if Guess.any():
        NPoints = len(Guess)
        DeltaX = (XN-X0)/NPoints
        U = Guess
    else:
        #NPoints = math.ceil((abs(XN-X0))/MaxStepSize)
        DeltaX = (XN-X0)/NPoints
        U = np.zeros(NPoints)
        
    #TrueStepSize = (XN-X0)/NPoints
    
    if Neuman and Robin:
        print("Error: Please specify only one type of boundary conditions")
        return 0
    elif Robin:
        delta = Bounds[1,1]
        UN = 2*delta*DeltaX
        gamma = Bounds[2,1]
        EndVector = (2,-1*(1+gamma*DeltaX))
    elif Neuman:
        delta = Bounds[1,1]
        UN = 2*delta*DeltaX
        EndVector = (2,-2)
    else:
        NPoints=NPoints-1
        U= U[:-1]
        UN = Bounds[1,1]
        EndVector = (1,-2)    
    
    XSpace = np.linspace(X0,XN,num=NPoints)
    #print(XSpace[2]-XSpace[1]+TrueStepSize)
    CoeffMatrix = np.zeros([NPoints,NPoints])
    CoeffMatrix[0,(0,1)] = -2,1
    CoeffVector = np.array([-1,0,1])
    for i in range(1,len(CoeffMatrix)-1):
        CoeffMatrix[i,i+CoeffVector] = (1,-2,1)
    CoeffMatrix[-1,(-2,-1)] = EndVector
    print(CoeffMatrix)
    
    #U = np.zeros(NPoints)
    U[0] = U0
    U[-1] = UN
    
    
        
    
    def EqnToSolve(SolnArray):
        return DiffusionConstant*(CoeffMatrix.dot(SolnArray)+U)/(DeltaX**2)+Reaction(XSpace,SolnArray)
        
    Ans = scipy.optimize.root(EqnToSolve,U)
    
    #Ans = Ans
    #print(Ans[0],Ans[1],Ans[0]-U[0],Ans[1]-Ans[0])
    #print(Ans[-1]+(Ans[1]-Ans[0]))
    if Ans.success:
        return XSpace, Ans.x
    else:
        print(Ans.message)
        return 0

def ExpectedSoln(BoundaryCond,Xarray,DiffusionConst=1):
    X0 = BoundaryCond[0,0]
    XN = BoundaryCond[0,1]
    U0 = BoundaryCond[1,0]
    UN = BoundaryCond[1,1]
    Soln = np.add((-1/(2*DiffusionConst))*np.multiply((Xarray-X0),(Xarray-XN)),(UN-U0)/(XN-X0)*(Xarray-X0)+U0)
    return Soln
#%%
if __name__ == "__main__":
    X0 = 0
    U0 = 0
    XN = 1
    UN = 0
    
    delta =1
    gamma =1
    
    DlechtBCs = np.array([[X0,XN],[U0,UN]]) 
    NeumanBCs = np.array([[X0,XN],[U0,delta]]) #For du/dx|Xn = delta
    RobinBCS = np.array([[X0,XN],[U0,delta],[0,gamma]]) #For du/dx|Xn = delta-gamma*u(Xn)
    print(RobinBCS)
    #%%
    #_,GuessForBratu = FiniteSolvePoisson(Bounds=BoundaryCond,Reaction=SimpleSourceTerm)
    XForGuess = np.linspace(X0,XN,num=40)
    GuessForBratu =ExpectedSoln(DlechtBCs,XForGuess)
    
    X,U = FiniteSolvePoisson(Bounds=DlechtBCs,Reaction=BratuTerm,Guess=GuessForBratu)
    

    plt.figure()
    plt.plot(X,U,label="My solution")
    ShowSolution = 0
    if ShowSolution:
        RealU = ExpectedSoln(BoundaryCond,X)
        plt.plot(X,RealU)
    plt.legend()
    plt.show()