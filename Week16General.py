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



def SingleShot(EqnToSolve,ShootArray):
   StartConditions = ShootArray[0:2]
   Period = ShootArray[-1]
   x,t = ODESolver.Solve_to(EqnToSolve,StartConditions,[0,Period],MinStep,ODESolver.RungeKutta4)
   G = StartConditions-x[:,-1] #Difference between output and input
   #print(G)
   G = np.append(G,EqnToSolve(x[:,-1],0)[0]) #Adds the 0 derrivative phase condition
   return G

ShootGuess = np.array([0.3,0.3,32])
def Shooting(EqnToSolve,ShootArray):
   CycleVector = scipy.optimize.root(lambda ShootArray: SingleShot(EqnToSolve,ShootArray),ShootArray)
   return CycleVector.x #find solution for a given array, should generalise with some modification

Test = SingleShot(predator_prey,ShootGuess)
print(Test,"worked")

FoundX,FoundY,FoundPeriod = Shooting(predator_prey,ShootGuess)
#print(Solution)
#FoundX,FoundY,FoundTime
FoundICs = [FoundX,FoundY]
FoundTime = [0, FoundPeriod]
x,t = ODESolver.Solve_to(predator_prey,FoundICs,FoundTime,MinStep,ODESolver.RungeKutta4)
plt.plot(x[0,:],x[1,:]) 
#plt.plot(t,x[1,:]) #Shows that period is roughly 32 for b= 0.1 and a 
plt.show()


def PeriodShooting(EqnToSolve,StartConditions,Period,Error,StepSize):
   OldG = np.linalg.norm(SingleShot(EqnToSolve, StartConditions,Period))
   if np.linalg.norm(OldG) > Error:
      NewG = np.linalg.norm(SingleShot(EqnToSolve, StartConditions,(Period+StepSize)))
      Period = Period+StepSize*(OldG/(NewG-OldG))
      return PeriodShooting(EqnToSolve,StartConditions,Period,Error,StepSize)
   else:
      return Period

#T = PeriodShooting(predator_prey,InitCon,[0,32],0.2)


#solution = scipy.optimize.fsolve(SingleShot)

### Any working solution would be using some weird method i dont at all undestand
### What does a solution mean




""" #Visualising the vector field
y0 = np.array([0.1,0.1])
t = np.linspace(0, 1, 20)
x,y = np.meshgrid(np.linspace(0.01,0.99,10),np.linspace(0.01,0.99,10))
EX , EY = predator_prey([x,y],0)
plt.streamplot(x,y,EX,EY, density=1.4, linewidth=None, color='#A23BEC')
plt.grid()
plt.show()
"""


#for i in range(1,10):
    #solution = scipy.integrate.odeint(predator_prey,y0*i,t)
    #plt.plot(i*y0[0],i*y0[1],'k*')
    #plt.plot(solution[:,0], solution[:,1])

   

#plt.plot(solution[:,0], solution[:,1]) #label="x")
#plt.xlabel("x")
#plt.plot(t, solution[:,1], label="y")
#plt.legend()
#plt.show()