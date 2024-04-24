# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 15:30:33 2024

@author: farka
"""

from CommonModules import *
import Week16General
import Week17Functions
import ODESolver
def SimpleSolveWrapper(Func,t,Param,SolnEstimate,SolverStepSize=0.001):
    WrappedFunc = lambda x: Func(x,t,Param)
    return scipy.optimize.root(WrappedFunc,SolnEstimate).x

def ShootingSolveWrapper(Func,t,Param,SolnEstimate,SolverStepSize=0.001):
    WrappedFunc = lambda u,t :Func(u, t,Param)
    x,period = Week16General.Shooting(WrappedFunc, SolnEstimate[:-1], SolnEstimate[-1],StepSize=SolverStepSize)
    return np.hstack([x,period])
    


def NaturalContinuation(Func,X0,Param0,WithShooting=0,ParamNSteps=10,ParamStepSize=0.1,SolverStepSize=0.1):
    """Natural Continuation

    Args:
        Func (function): function for continuation
        X0 (array): initial conditions
        Param0 (scalar): starting parameter value
        WithShooting (boolean): option to use shooting methods, default is off
        ParamNSteps (int, optional): Number of continuation steps. Defaults to 10.
        ParamStepSize (float, optional): Size of each step, can be positive or negative. Defaults to 0.1.
        SolverStepSize (float, optional): Step size used by the ODE solver. Defaults to 0.1.

    Returns:
        SolutionSpace(2d-array) : array of solutions for each parameter value
        ParamSpace(1d-array) : corresponding array of parameter values
    """
    ParameterSpace =np.linspace(Param0,Param0+ParamStepSize*ParamNSteps,ParamNSteps)
    SolutionSpace = np.zeros([ParamNSteps,(np.size(X0))])
    if WithShooting:
        WrapperToUse = ShootingSolveWrapper
    else:
        WrapperToUse = SimpleSolveWrapper
    SolutionSpace[-1,:] = X0 #Uses last part of solution as temporary storage for initial approximation
    for i in range(0,ParamNSteps):
        try:
            SolutionSpace[i,:] = WrapperToUse(Func,1,ParameterSpace[i],SolutionSpace[i-1,:],SolverStepSize=SolverStepSize)
        except:
            if i:
                print(f"Failed to go past parameter= {ParameterSpace[i]}: Returned all successful attempts")
                return SolutionSpace[:i-1,:] ,ParameterSpace[:i-1]
            else:
                print(f"Continuation Failed")
                return 0,0
    return SolutionSpace ,ParameterSpace




def ShootingArcLengthCont(Func,X0,ParamBounds,ContinuationMaxSteps,
                        ParamStepSize=0.1,Solver=ODESolver.RungeKutta4,
                        SolverStepSize=0.1
                        ):
    """_summary_

    Args:
        Func (_type_): _description_
        X0 (_type_): array of [x1...xn,t]
        ParamBounds (tuple): range of parameter to check across (A,B)
        ContinuationMaxSteps (_type_): max number of iterations
        ParamStepSize (float, optional): First step in parameter space, Defaults to 0.1.
        Solver (_type_, optional): integrator to use, Defaults to ODESolver.RungeKutta4.
        SolverStepSize (float, optional): stepsize of numerical integrator, Defaults to 0.1.

    Returns:
        SolnSpace: 2d array of solutions
        ParamSpace: 1d array of corresponding parameters
    """
    P0,PN = ParamBounds[0],ParamBounds[1]
    Direction = np.sign(PN-P0)
    #ParamArray = np.linspace(P0,PN,num=MaxSteps)
    ParamSpace = np.empty(ContinuationMaxSteps)+np.nan
    SolnSpace = np.empty((ContinuationMaxSteps,np.size(X0)))+np.nan
    
    ParamSpace[0],ParamSpace[1] = ParamBounds[0],ParamBounds[0]+Direction*ParamStepSize
    SolnSpace[0,:] = ShootingSolveWrapper(Func,1,ParamSpace[0],X0)
    SolnSpace[1,:] = ShootingSolveWrapper(Func,1,ParamSpace[1],SolnSpace[0,:])
    
    def ArcShootRootFind(SolnAndParam):
        
        SolnToInvestigate,ParamToInvestigate = SolnAndParam[:-1],SolnAndParam[-1]
        FuncAtParam = lambda x,t : Func(x,t,ParamToInvestigate)
        
        ShootingLHS = Week16General.SingleShot(SolnToInvestigate,SolverStepSize,Solver,FuncAtParam)
        ArcCondition = (np.dot(ParamSecant,(ParamToInvestigate-ParamEstimate))+
                                np.dot(SolnSecant,(SolnToInvestigate-SolnEstimate)))
        ArcLengthRHS = np.append(ShootingLHS,ArcCondition)
        
        return ArcLengthRHS
    
    for i in range(2,ContinuationMaxSteps):
            #SolnSecant =Soln[i-1]-SolnOlder
            #ParamSecant =ParamOld-ParamOlder # this could be done in the loop, would make it easier
            ParamSecant = ParamSpace[i-1]-ParamSpace[i-2]
            SolnSecant = SolnSpace[i-1,:]-SolnSpace[i-2,:]
            SolnEstimate = SolnSpace[i-1,:]+SolnSecant
            ParamEstimate = ParamSpace[i-1]+ParamSecant
            #print(SolnEstimate,"soln")
            #print(ParamEstimate,"Param")
            try:
                Result = scipy.optimize.root(ArcShootRootFind,np.hstack((SolnEstimate,ParamEstimate))).x
            except:
                SolnSpace = np.squeeze([i for i in SolnSpace if not np.isnan(i).all()])
                ParamSpace = [i for i in ParamSpace if not np.isnan(i)]
                print(f"Continuation Failed after {i} steps")
                return SolnSpace,ParamSpace
                #break
            #NewSoln,NewParam = Result[0],Result[1]
            SolnSpace[i,:] = Result[:-1]
            ParamSpace[i] = Result[-1]
            

    SolnSpace = np.squeeze([i for i in SolnSpace if not np.isnan(i).all()])
    ParamSpace = [i for i in ParamSpace if not np.isnan(i)]    
    return SolnSpace,ParamSpace
    #Would need to do natural for first 2
    
    
    
    
    



        
def Main():
    #Answer,Space = NaturalParameterContintuation(Week17Functions.Cubic, 1.1, -2, 40) 
    Answer,Space = ShootingNumericalContinuation(Week17Functions.ModifiedBetaFormHopf,[2.3,0,30],2,ParamNSteps=2,ParamStepSize=-0.1 )
    plt.plot(Space,Answer[:,0])
    plt.show()
    return       



if __name__ == "__main__":
    Main()
    

