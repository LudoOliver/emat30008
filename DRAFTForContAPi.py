"""Testing Generalisation of Continuation"""
from CommonModules import *
from Week17Functions import *
import ODESolver
import Week17Continuation
def Continuation(Func,X0,Param0,Shooting=False,ArcLength=False,ParamStepSize=0.1,ParamEndPoint=None,ParamNSteps=None):
    """_summary_

    Args:
        Eqn (function): function of the form f(x,t,Parameter)
        X0 (array): initial conditions
        Param0 (float): parameter start value
        Shooting (int, optional): Shooting on/off flag. Defaults to 0.
        ArcLength (int, optional): ArcLength on/off. Defaults to 0.
    """
    if ParamEndPoint==None and ParamNSteps==None:
        print("Please choose either an end point or number of iterations ")
    if not Shooting and not ArcLength:
        return Week17Continuation.NaturalParameterContintuation(Func,X0,Param0,ParamNSteps=ParamNSteps,ParamStepSize=ParamStepSize)
    elif Shooting and not ArcLength:
        return Week17Continuation.ShootingNumericalContinuation(Func,X0,Param0,ParamNSteps=ParamNSteps,ParamStepSize=ParamStepSize)
    elif Shooting and ArcLength:
        
    