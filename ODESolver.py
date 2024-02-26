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

#def Solve_to(t0,t1,x0,DeltaTMax,SolverToUse,FuncToSolve):
def Solve_to(FuncToSolve,x0,tspan,DeltaTMax,SolverToUse,):
    #print(tspan)
    t0,t1 = tspan
    NSteps = math.ceil((t1-t0)/DeltaTMax) #Whole Number of step size
    StepSize = (t1-t0)/NSteps #Creates Step size less than maximum
    TArray = np.linspace(t0,t1,NSteps)
    XArray = np.zeros([np.size(x0),NSteps])
    XArray[:,0] = x0
    for i in range(1,NSteps):
        XArray[:,i] = SolverToUse(FuncToSolve,XArray[:,i-1],TArray[i-1],StepSize)
    return XArray, TArray

#%% Specific
"""
def ODEFunc(t,x):
    return x

def VectorODe(t,x):
    x1 = x[1]
    x2 = -x[0]
    return np.array([x1,x2])
x0 = [0,1]
t0 = 0
t1 = 50
Time = [0,50]
MinStep = 0.1

t,x = Solve_to(VectorODe,x0,Time,MinStep,EulerStep,)
plt.subplot(2,1,1)
plt.plot(t,x[0,:])
plt.title("x vs t")
plt.subplot(2,1,2)
plt.plot(x[1,:],x[0,:])
plt.title("x vs x'")
plt.tight_layout()
plt.suptitle("Eulers method solutions", fontsize=16)
plt.subplots_adjust(top=0.85,bottom=0.15)
plt.figtext(0.1,0.05,f"For x''=-x with x(0)={x0[0]} and x'(0)={x0[1]}  evaluated to x({t1})with stepsize {MinStep}",wrap=True)

plt.show()

#plt.plot(x,y[],label=f"RungeKutta4")
#t,x = Solve_to(0,1,1,0.001,EulerStep,ODEFunc)
#plt.plot(x,y,label=f"Eulers method")
#plt.legend()
#plt.show()

"""