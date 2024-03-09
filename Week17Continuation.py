# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 15:30:33 2024

@author: farka
"""

from CommonModules import *
import Week16General
import Week17Functions

def SimpleSolveWrapper(Func,t,Param,SolnEstimate):
    WrappedFunc = lambda x: Func(x,t,Param)
    return scipy.optimize.root(WrappedFunc,SolnEstimate).x

def ShootingSolveWrapper(Func,t,Param,SolnEstimate):
    WrappedFunc = lambda u,t :Func(u, t,Param)
    #print(SolnEstimate[:-1])
    #print(SolnEstimate[-1])
    print(Param)
    #print(SolnEstimate)
    x,period = Week16General.Shooting(WrappedFunc, SolnEstimate[:-1], SolnEstimate[-1])#,StepSize=0.0001)
    return np.hstack([x,period])
    
def NaturalParameterContintuation(Func,X0,Param0,ParamNSteps,ParamStep=0.1):#,discretisation=lambda x:x):
    ParameterSpace =np.linspace(Param0,Param0+ParamStep*ParamNSteps,ParamNSteps)
    SolutionSpace = np.zeros([ParamNSteps,(np.size(X0))])
    SolutionSpace[0,:] = scipy.optimize.root(lambda x: Func(x,1,Param0),X0).x
    for i in range(1,ParamNSteps):
        #CurrentParam = ParameterSpace[i]
        SolutionSpace[i,:] = SimpleSolveWrapper(Func,1,ParameterSpace[i],SolutionSpace[i-1,:])
    return SolutionSpace, ParameterSpace 

    
def ShootingNumericalContinuation(Func,X0,Param0,ParamNSteps=10,ParamStepSize=0.1):
    ParameterSpace =np.linspace(Param0,Param0+ParamStepSize*ParamNSteps,ParamNSteps)
    SolutionSpace = np.zeros([ParamNSteps,(np.size(X0))])
    Initial = lambda u,t :Func(u, t,Param0)
    SolutionSpace[0,:-1],SolutionSpace[0,-1] = Week16General.Shooting(Initial, X0[0:-1], X0[-1]) #could be unpack issues
    for i in range(1,ParamNSteps):
        #Current = lambda u,t :Func(u, t,ParameterSpace[i])
        #SolutionSpace[i,:-1],SolutionSpace[i,-1] = Week16General.Shooting(Current, SolutionSpace[i-1,:-1], SolutionSpace[i-1,-1])
        SolutionSpace[i,:] = ShootingSolveWrapper(Func,1,ParameterSpace[i],SolutionSpace[i-1,:])
        
    
    return SolutionSpace ,ParameterSpace

def OneArcLengthCont(func,SolnArray,ParamArray):
    ParamSecant = ParamArray[-1]-ParamArray[-2]
    SolnArray = SolnArray[:,-1]-ParamArray[:,-2]
    



#def ArcLengthContinuation
    
        
def Main():
    #Answer = NaturalParameterContintuation(Week17Functions.Cubic, 1.1, -2, 40) 
    Answer,Space = ShootingNumericalContinuation(Week17Functions.ModifiedBetaFormHopf,[2.3,0,30],2,ParamNSteps=2,ParamStepSize=-0.1 )
    plt.plot(Space,Answer[:,0])
    plt.show()
    return       



if __name__ == "__main__":
    Main()