"""1st order solver"""
import math
import numpy as np
import matplotlib.pyplot as plt
#%% General

def EulerStep(Ftx, Xn, Tn,StepSize=0.001):
    return np.add(Xn,StepSize*Ftx(Xn,Tn))

def RungeKutta4(Ftx,Xn,Tn,StepSize=0.001):
    k1 = Ftx(Xn,Tn)
    k2 = Ftx(Xn+StepSize*k1/2,Tn+(StepSize/2))
    k3 = Ftx(Xn+StepSize*k2/2,Tn+StepSize/2)
    k4 = Ftx(Xn+StepSize*k3,Tn+StepSize)
    return np.add(Xn,(StepSize/6)*(k1+2*k2+2*k3+k4))


def Solve_to(FuncToSolve,x0,tspan,DeltaTMax=0.0001,SolverToUse=RungeKutta4):#(Ftx, Xn, Tn)):
    """ Solves a given ODE
    Arguements:
        FuncToSolve (function(x,t)): rhs of a ode
        x0 (array): array of initial conditions
        tspan (tuple): start and end of timespan to integrate across
        DeltaTMax (float): maximum timestep to use, default = 0.001
        SolverToUse (function): integrator to use, default is RungeKutta4
        

    Returns:
        X (2d array): solutions of form x(t) to the ODE
        T (1d array): corresponding time values
    """
    t0,t1 = tspan  # unpacking the time values
    NSteps = math.ceil((t1-t0)/DeltaTMax) #Whole Number of step size
    StepSize = (t1-t0)/NSteps #Creates Step size less than maximum
    TArray = np.linspace(t0,t1,NSteps)
    XArray = np.zeros([np.size(x0),NSteps])
    XArray[:,0] = x0
    for i in range(1,NSteps):
        XArray[:,i] = SolverToUse(FuncToSolve,XArray[:,i-1],TArray[i-1],StepSize)
    return XArray, TArray
