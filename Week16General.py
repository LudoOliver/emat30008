
import scipy
import math
import numpy as np
import matplotlib.pyplot as plt
import ODESolver
from scipy.optimize import root

def predator_prey(x,t):
   #print(x)
   dx = (x[0])*(1-(x[0])) - ( (a*(x[0])*(x[1])) / (d+(x[0])) )
   dy = ( b*x[1] )*(1 - ( x[1] / x[0]))
   return np.array([dx , dy])
#print((6)(1))




def SingleShot(ShootArray,StepSize,Solver,EqnToSolve):
    """The root finding problem

    Args:
        ShootArray ([Array],Scalar): X0, Time period
        StepSize (float): Solver step size
        Solver (function): ode solver to use
        EqnToSolve (function): chosen numerical integrator

    Returns:
        Numpy array of X(0)-X(T) stacked upon dx_1/dt
    """
    StartConditions = ShootArray[0:-1]
    Period = ShootArray[-1]
    x,t = ODESolver.Solve_to(EqnToSolve,StartConditions,[0,Period],StepSize,Solver)
    G = StartConditions-x[:,-1] #Difference between output and input
    G = np.append(G,EqnToSolve(x[:,-1],0)[0]) #Adds the 0 derrivative phase condition
    return G

def Shooting(EqnToSolve,X0,T0,StepSize=0.001,Solver=ODESolver.RungeKutta4):
    """
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
        print("Error: Dimension of X0 is too large")
        return 0
    try:    
        #CycleVector = scipy.optimize.root(lambda ShootArray: SingleShot(ShootArray,StepSize,Solver),ShootArray)
        #NewVector = root(SingleShot,InitialGuess,args=SingleShotArgs).x
        SolnVec = root(NewSingleShot,InitialGuess,args=(StepSize,Solver,EqnToSolve)).x#,args=SingleShotArgs).x
    except:
        print("Shooting Error: Try Checking Your Initial Condition or Decreasing Your Stepsize")
        return 0
        #print(CycleVector.x)
    return SolnVec[:-1],SolnVec[-1] #



#%%
def main():
    a = 1
    d = 0.1
    b = 0.1
    InitCon = np.array([0.1,0.1])
    MinStep = 0.1
    Time = [0,100]
    ShootGuess = np.array([0.3,0.3,32])
    #ShootX = np.array([0.1,0.1])
    ShootX = np.array([1,3,4,5])
    ShootT=32
    FoundX,FoundPeriod = Shooting(predator_prey,ShootX,ShootT)
# print(Solution)
# FoundX,FoundY,FoundTime
    FoundICs = FoundX
    FoundTime = [0, FoundPeriod]
    x,t = ODESolver.Solve_to(predator_prey,FoundICs,FoundTime,MinStep,ODESolver.RungeKutta4)
    plt.plot(x[0,:],x[1,:]) 
    #plt.plot(t,x[1,:]) #Shows that period is roughly 32 for b= 0.1 and a 
    #plt.plot(t,x[0,:])
    plt.show()
if __name__ =="__main__":
    a = 1
    d = 0.1
    b = 0.1
    InitCon = np.array([0.1,0.1])
    MinStep = 0.1
    Time = [0,100]
    ShootGuess = np.array([0.3,0.3,32])
    ShootX = np.array([0.1,0.1])
    main()

""" Junk code, kept for testing /if i need to go back quickly
    def SingleShot(ShootArray,StepSize=StepSize,Solver=Solver):
        StartConditions = ShootArray[0:-1]
        Period = ShootArray[-1]
        x,t = ODESolver.Solve_to(EqnToSolve,StartConditions,[0,Period],StepSize,Solver)
        G = StartConditions-x[:,-1] #Difference between output and input
        G = np.append(G,EqnToSolve(x[:,-1],0)[0]) #Adds the 0 derrivative phase condition
        return G
    SingleShotArgs = (StepSize,Solver)
def SingleShot(EqnToSolve,ShootArray,StepSize,Solver):
   StartConditions = ShootArray[0:-1]
   Period = ShootArray[-1]
   x,t = ODESolver.Solve_to(EqnToSolve,StartConditions,[0,Period],StepSize,Solver)
   G = StartConditions-x[:,-1] #Difference between output and input
   #print(G)
   G = np.append(G,EqnToSolve(x[:,-1],0)[0]) #Adds the 0 derrivative phase condition
   return G
"""