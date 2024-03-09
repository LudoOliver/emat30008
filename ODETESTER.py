"""Investigating ODEs"""
from CommonModules import *
from Week17Functions import *
import ODESolver
Param0 = 2
FixBeta = lambda u,t :ModifiedBetaFormHopf(u, t,Param0)
x,t = ODESolver.Solve_to(FixBeta,[1.1,1.1],[0,30])

import numpy as np
import matplotlib.pyplot as plt
#%matplotlib inline

x,y = np.meshgrid(np.linspace(-3,3,100),np.linspace(-3,3,100))

u = -y/np.sqrt(x**2 + y**2)
v = x/np.sqrt(x**2 + y**2)
[u,v] =  FixBeta([x,y],1)

plt.quiver(x,y,u,v)
plt.show()
#%%

