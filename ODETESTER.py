"""Investigating ODEs"""
from CommonModules import *
from Week17Functions import *

import ODESolver
Param0 = 1
FixBeta = lambda u,t :BetaFormHopf(u, t,Param0)
x,t = ODESolver.Solve_to(FixBeta,[1,0.1],[0,6])
plt.plot(t,x[0,:])
plt.plot(t,x[1,:])
plt.title(f"TimeBehaviour for beta = {Param0}")
plt.show()
plt.plot(x[0,:],x[1,:])
plt.title("phase behaviour in x,y space")
import numpy as np
import matplotlib.pyplot as plt
#%matplotlib inline
plt.show()
x,y = np.meshgrid(np.linspace(-3,3,100),np.linspace(-3,3,100))

u = -y/np.sqrt(x**2 + y**2)
v = x/np.sqrt(x**2 + y**2)
[u,v] =  FixBeta([x,y],1)

plt.quiver(x,y,u,v)
plt.show()
#%%
print(1+[1,2])


# %%
