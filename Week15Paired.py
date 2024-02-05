# Problem 1

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
MinStep = 0.01
Time = [0,100]

x,t = ODESolver.Solve_to(predator_prey,InitCon,Time,MinStep,ODESolver.RungeKutta4)
#print(x)
#print(t)
plt.plot(x[0,:],x[1,:])
plt.show()
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