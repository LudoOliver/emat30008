
import scipy
import math
import numpy as np
import matplotlib.pyplot as plt
import ODESolver
from scipy.optimize import root

def SingleShot(ShootArray,StepSize,Solver,EqnToSolve):
    """Define the root finding problem for shooting methods, not intended as a user level function

    Args:
        ShootArray ([Array],Scalar): X0, Time period
        StepSize (float): Solver step size
        Solver (function): ODE solver to use
        EqnToSolve (function): chosen numerical integrator

    Returns:
        Root finding problem for shooting, with zero derivative phase condition
    """
    StartConditions = ShootArray[0:-1]
    Period = ShootArray[-1]
    x,t = ODESolver.Solve_to(EqnToSolve,StartConditions,[0,Period],StepSize,Solver)
    G = StartConditions-x[:,-1] #Difference between output and input
    G = np.append(G,EqnToSolve(x[:,-1],0)[0]) #Adds the 0 derrivative phase condition
    return G

def Shooting(EqnToSolve,X0,T0,StepSize=0.001,Solver=ODESolver.RungeKutta4):
    """
    Preforms numerical shooting on a given ODE
    Input:
        EqnToSolve : Function of(x,t)
            -the ODE to integrate
        X0 : Array
            -the initial conditions
            -should be as close to the suspected limit cycle as possible
        T : Scalar
            - an approxmiation of the cycles period
        Solver: Function
            -Numerical Integrator to use  
            -Must only take arguements (function, x, t, StepSize)
            -default is RungeKutta4
        StepSize:
            -Step size used in for numerical integration
            -default = 0.001
    Output:
        X: Array matching dimensions X0
            -the location of the limit cycle in x space
        T: Scalar    
            -the time period of the limit   
    """
    ShootArray = np.append(X0,T0)
    InitialGuess = np.append(X0,T0)
    try:
        OutputDim= EqnToSolve(X0,1)
    except:
        print("Error: Invalid X0")
        return 0
    if np.size(OutputDim)!=np.size(X0):
        print("Error: Incorrect dimensions of X0")
        return 0
    try:    
        SolnVec = root(SingleShot,InitialGuess,args=(StepSize,Solver,EqnToSolve)).x
    except:
        print("Shooting Error: Try Checking Your Initial Condition or Decreasing Your Stepsize")
        return 0
        #print(CycleVector.x)
    return SolnVec[:-1],SolnVec[-1] #


