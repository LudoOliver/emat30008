# Problem 1

### Main problems - How do i work through this ive got 3 unkowns im trying to find roots of right?
### unlesss i x corrrealte to find period wtf do i do
### do i just care about finding a single co-ordinate?
### would it be easisest to just iterate for ages and find repeated section in

import scipy
import math
import numpy as np
import matplotlib.pyplot as plt
import ODESolver
a = 1
d = 0.1
b = 0.1

def predator_prey(x,t):
   #print(x)
   dx = (x[0])*(1-(x[0])) - ( (a*(x[0])*(x[1])) / (d+(x[0])) )
   dy = ( b*x[1] )*(1 - ( x[1] / x[0]))
   return np.array([dx , dy])
#print((6)(1))

InitCon = np.array([0.1,0.1])
MinStep = 0.1
Time = [0,100]



def SingleShot(EqnToSolve,ShootArray,StepSize,Solver):
   StartConditions = ShootArray[0:-1]
   Period = ShootArray[-1]
   x,t = ODESolver.Solve_to(EqnToSolve,StartConditions,[0,Period],StepSize,Solver)
   G = StartConditions-x[:,-1] #Difference between output and input
   #print(G)
   G = np.append(G,EqnToSolve(x[:,-1],0)[0]) #Adds the 0 derrivative phase condition
   return G

ShootGuess = np.array([0.3,0.3,32])
ShootX = np.array([0.1,0.1])
ShootT=32

def Shooting(EqnToSolve,X0,T0,StepSize=0.001,Solver=ODESolver.RungeKutta4):
   """
    Input:
        EqnToSolve : Function of(x,t)
            -the ODE to integrate
        X0 : Array
            -the initial conditions
            -should be in close to the suspcted limit cycle
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
   try:
       CycleVector = scipy.optimize.root(lambda ShootArray: SingleShot(EqnToSolve,ShootArray,StepSize,Solver),ShootArray)
       #print("Success")
   except:
       print("Error: Try Checking Your Initial Condition or Decreasing Your Stepsize")
       return 
    #print(CycleVector.x)
   return CycleVector.x[:-1],CycleVector.x[-1] #find solution for a given array, should generalise with some modification

#Test = SingleShot(predator_prey,ShootGuess)
#print(Test,"worked")



#%%
def main():
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
    main()

