# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 15:30:33 2024

@author: farka
"""

from CommonModules import *
import Week16General
import Week17Functions
    
def NaturalParameterContintuation(Func,X0,Param0,ParamNSteps,ParamStep=0.1):#,discretisation=lambda x:x):
    ParameterSpace =np.linspace(Param0,Param0+ParamStep*ParamNSteps,ParamNSteps)
    SolutionSpace = np.zeros([ParamNSteps,(np.size(X0))])
    SolutionSpace[0,:] = scipy.optimize.root(lambda x: Func(x,1,Param0),X0).x
    for i in range(1,ParamNSteps):
        CurrentParam = ParameterSpace[i]
        #print(CurrentParam)
        Guess = scipy.optimize.root(lambda x: Func(x,1,CurrentParam),SolutionSpace[i-1,:])
        #print(Guess.success)
        SolutionSpace[i,:] = Guess.x
    return SolutionSpace, ParameterSpace 

# def ShootingParameterContintuation(Func,X0,Param0,ParamNSteps,ParamStep=0.1):#,discretisation=lambda x:x):
#     ParameterSpace =np.linspace(Param0,Param0+ParamStep*ParamNSteps,ParamNSteps)
#     SolutionSpace = np.zeros([ParamNSteps,(np.size(X0))])
#     SolutionSpace[0,:] = scipy.optimize.root(lambda x: Func(x,1,Param0),X0).x
#     for i in range(1,ParamNSteps):
#         CurrentParam = ParameterSpace[i]
#         #print(CurrentParam)
#         #Guess = scipy.optimize.root(lambda x: Func(x,1,CurrentParam),SolutionSpace[i-1,:])
#         Guess = Week16General.Shooting(Func, SolutionSpace[i-1,:-1], SolutionSpace[i-1,-1])
#         #print(Guess.success)
#         #SolutionSpace[i,:] = Guess.x
#     return SolutionSpace, ParameterSpace 
#def ContinuationWrapper(ode,method="natural"):
    
def ShootingNumericalContinuation(func,X0,Param0,ParamNSteps=10,ParamStepSize=0.1):
    ParameterSpace =np.linspace(Param0,Param0+ParamStepSize*ParamNSteps,ParamNSteps)
    SolutionSpace = np.zeros([ParamNSteps,(np.size(X0))])
    Initial = lambda u,t :func(u, t,Param0)
    SolutionSpace[0,:-1],SolutionSpace[0,-1] = Week16General.Shooting(Initial, X0[0:-1], X0[-1]) #could be unpack issues
    for i in range(1,ParamNSteps):
        Current = lambda u,t :func(u, t,ParameterSpace[i])
        SolutionSpace[i,:-1],SolutionSpace[i,-1] = Week16General.Shooting(Current, SolutionSpace[i-1,:-1], SolutionSpace[i-1,-1])
        
    return SolutionSpace ,ParameterSpace

#def ArcLengthContinuation
    
        
def Main():
    Answer = NaturalParameterContintuation(Week17Functions.Cubic, 1.1, -2, 40) 
    return Answer       

if __name__ == "__main__":
    Answer,Space = Main()
    for i in range(len(Answer)):
        plt.scatter(Space[i],Week17Functions.Cubic(Answer[i], 0,Space[i]))
