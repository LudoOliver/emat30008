"""1st order solver"""
import math
import numpy as np
import matplotlib.pyplot as plt
#%% General
def EulerStep(Tn,Xn,StepSize,Ftx):
    return Xn+StepSize*Ftx(Xn,Tn)

def Solve_to(t1,t2,x1,DeltaTMax,SolverToUse,FuncToSolve):
    NSteps = math.ceil((t2-t1)/DeltaTMax) #Whole Number of step size
    StepSize = (t2-t1)/NSteps #Creates Step size less than maximum
    TArray = np.linspace(t1,t2,NSteps)
    XArray = np.zeros(NSteps)
    XArray[0] = x1
    for i in range(1,NSteps):
        XArray[i] = SolverToUse(TArray[i-1],XArray[i-1],StepSize,FuncToSolve)
    return TArray,XArray

#%% Specific
def ODEFunc(x,t):
    return x
for i in [0.1,0.01,0.001,0.0001]:
    x,y = Solve_to(0,1,1,i,ODEFunc)
    plt.plot(x,y,label=f"MaxStep={i}")
plt.legend()
plt.show()