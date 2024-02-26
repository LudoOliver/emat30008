# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 15:30:33 2024

@author: farka
"""

import CommonModules
import Week16General
import Week17Functions
    
def NaturalParameterContintuation(Func,X0,Param0,ParamNSteps,ParamStep=0.1):#,discretisation=lambda x:x):
    ParameterSpace =np.linspace(Param0,Param0+ParamStep*ParamNSteps,ParamNSteps)
    SolutionSpace = np.zeros((size(X0),ParamNSteps))
    SolutionSpace[0,:] = X0
    for i in range(1,ParamNSteps):
        CurrentParam = ParameterSpace[i]
        SolutionSpace[i,:] = scipy.optimize.root(lambda x: Func(x,CurrentParam),SolutionSpace[i:-1])
        
    
    
    