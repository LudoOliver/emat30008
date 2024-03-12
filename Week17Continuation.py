# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 15:30:33 2024

@author: farka
"""

from CommonModules import *
import Week16General
import Week17Functions
import ODESolver
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
    
def NaturalParameterContintuation(Func,X0,Param0,ParamNSteps,ParamStepSize=0.1):#,discretisation=lambda x:x):
    ParameterSpace =np.linspace(Param0,Param0+ParamStepSize*ParamNSteps,ParamNSteps)
    SolutionSpace = np.zeros([ParamNSteps,(np.size(X0))])
    SolutionSpace[0,:] = scipy.optimize.root(lambda x: Func(x,1,Param0),X0).x
    for i in range(1,ParamNSteps):
        #CurrentParam = ParameterSpace[i]
        SolutionSpace[i,:] = SimpleSolveWrapper(Func,1,ParameterSpace[i],SolutionSpace[i-1,:])
    return SolutionSpace, ParameterSpace 

    
def ShootingNumericalContinuation(Func,X0,Param0,ParamNSteps=10,ParamStepSize=0.1,SolverStepSize=0.1):
    #0.1 Is recomended as it is incredibly slow for Modified Btea From
    ParameterSpace =np.linspace(Param0,Param0+ParamStepSize*ParamNSteps,ParamNSteps)
    SolutionSpace = np.zeros([ParamNSteps,(np.size(X0))])
    Initial = lambda u,t :Func(u, t,Param0)
    SolutionSpace[0,:-1],SolutionSpace[0,-1] = Week16General.Shooting(Initial, X0[0:-1], X0[-1],SolverStepSize=0.1) #could be unpack issues
    for i in range(1,ParamNSteps):
        #Current = lambda u,t :Func(u, t,ParameterSpace[i])
        #SolutionSpace[i,:-1],SolutionSpace[i,-1] = Week16General.Shooting(Current, SolutionSpace[i-1,:-1], SolutionSpace[i-1,-1])
        SolutionSpace[i,:] = ShootingSolveWrapper(Func,1,ParameterSpace[i],SolutionSpace[i-1,:])
        
    
    return SolutionSpace ,ParameterSpace

def NaturalContinuation(Func,X0,Param0,WithShooting=0,ParamNSteps=10,ParamStepSize=0.1,SolverStepSize=0.1):
    ParameterSpace =np.linspace(Param0,Param0+ParamStepSize*ParamNSteps,ParamNSteps)
    SolutionSpace = np.zeros([ParamNSteps,(np.size(X0))])
    if WithShooting:
        WrapperToUse = ShootingSolveWrapper
    else:
        WrapperToUse = SimpleSolveWrapper
    #Could be clever putting IC's to front?
    SolutionSpace[-1,:] = X0
    for i in range(0,ParamNSteps):
        SolutionSpace[i,:] = WrapperToUse(Func,1,ParameterSpace[i],SolutionSpace[i-1,:])
    return SolutionSpace ,ParameterSpace

#def ShootingArcLengthWrapper(func,X,Param)
# return func



def ShootingArcLengthCont(Func,X0,ParamBounds,ContinuationMaxSteps,
                        ParamStepSize=0.1,WithShooting=1,Solver=ODESolver.RungeKutta4,
                        SolverStepSize=0.1
                        ):
    
    P0,PN = ParamBounds[0],ParamBounds[1]
    Direction = np.sign(PN-P0)
    #ParamArray = np.linspace(P0,PN,num=MaxSteps)
    ParamSpace = np.empty(MaxSteps)+np.nan
    SolnSpace = np.empty(np.size(X0),MaxSteps)+np.nan
    
    ParamSpace[0],ParamSpace[1] = ParamBounds[0],ParamBounds[0]+Direction*ParamStepSize
    SolnSpace[0,:] = ShootingSolveWrapper(Func,1,ParamSpace[0],X0)
    SolnSpace[1,:] = ShootingSolveWrapper(Func,1,ParamSpace[1],SolnSpace[0,:])
    
    def OldArcShootRootFind(SolnAndParam):
        
        SolnToInvestigate,ParamToInvestigate = SolnAndParam[0],SolnAndParam[1]
        Period = SolnToInvestigate[-1]
        
        FuncAtParam = lambda x,t : Func(x,t,ParamToInvestigate)
        X0 = SolnEstimate[:-1]
        X,_ = ODESolver.Solve_to(Func,X0,[0,Period],DeltaTMax=0.1)
        EndX = X[-1,:]
        XCondition = X-X0
        PhaseCondition = FuncAtParam(X0)[0]
        ArcCondition = (np.dot(ParamSecant,(ParamToInvestigate-ParamEstimate))+
                                np.dot(SolnSecant,(SolnToInvestigate-SolnEstimate)))
        
        return np.array([XCondition,PhaseCondition,ArcCondition])
    
    def ArcShootRootFind(SolnAndParam):
        
        SolnToInvestigate,ParamToInvestigate = SolnAndParam[0],SolnAndParam[1]
        FuncAtParam = lambda x,t : Func(x,t,ParamToInvestigate)
        SingleShotArgs = (SolverStepSize,Solver,FuncAtParam)
        
        ShootingLHS = Week16General.SingleShot(SolnToInvestigate,SolverStepSize,Solver,FuncAtParam)
        ArcCondition = (np.dot(ParamSecant,(ParamToInvestigate-ParamEstimate))+
                                np.dot(SolnSecant,(SolnToInvestigate-SolnEstimate)))
        ArcLengthRHS = np.append(ShootingLHS,ArcCondition)
        
        return ArcLengthRHS
    
    for i in range(2,ContinuationMaxSteps):
        if min(ParamBounds) < ParamSpace[i-1] < max(ParamBounds):
            #SolnSecant =Soln[i-1]-SolnOlder
            #ParamSecant =ParamOld-ParamOlder # this could be done in the loop, would make it easier
            ParamSecant = ParamSpace[i-1]-ParamSpace[i-2]
            SolnSecant = SolnSpace[i-1,:]-SolnSpace[i-2,:]
            SolnEstimate = SolnSpace[i-1,:]+SolnSecant
            ParamEstimate = ParamSpace[i-1]+ParamSecant
            
            Result = scipy.optimize.root(ArcShootRootFind,[SolnEstimate,ParamEstimate])
            #NewSoln,NewParam = Result[0],Result[1]
            SolnSpace[i,:] = Result[0]
            ParamSpace[i] = Result[1]
            
        else:
            SolnSpace = np.squeeze([i for i in SolnSpace if not np.isnan(i).all()])
            ParamSpace = [i for i in ParamSpace if not np.isnan(i)]
            print("Bounds Exceeded")
            return SolnSpace,ParamSpace
    SolnSpace = np.squeeze([i for i in SolnSpace if not np.isnan(i).all()])
    ParamSpace = [i for i in ParamSpace if not np.isnan(i)]    
    return SolnSpace,ParamSpace
    #Would need to do natural for first 2
    
    
    
    
    
    
    
def ArcShootingSolveWrapper(Func,t,ParamArray,SolnArray):
    #Root finding problem
    
    
    
    
    
    
    WrappedFunc = lambda u,t :Func(u, t,Param)
    #print(SolnEstimate[:-1])
    #print(SolnEstimate[-1])
    print(Param)
    #print(SolnEstimate)
    x,period = Week16General.Shooting(WrappedFunc, SolnEstimate[:-1], SolnEstimate[-1])#,StepSize=0.0001)
    CurrentParam = #SOmething
    CurrentSoln = #Something
    #return np.hstack([x,period])    
    ParamSecant = ParamArray[-1]-ParamArray[-2]
    SolnSecant = SolnArray[:,-1]-SolnArray[:,-2]
    ParamCondition = np.dot(ParamSecant,(CurrentParam-ParamArray[-1]-ParamSecant))
    SolnCondition = np.dot(SolnSecant,(CurrentSoln-SolnArray[-1]-SolnSecant))
    ArcConditon = ParamCondition+SolnCondition                 
    x,period = Week16General.Shooting(WrappedFunc, SolnEstimate[:-1], SolnEstimate[-1])#,StepSize=0.0001)
    
    return np
    
    
    



#def ArcLengthContinuation
    
        
def Main():
    #Answer,Space = NaturalParameterContintuation(Week17Functions.Cubic, 1.1, -2, 40) 
    Answer,Space = ShootingNumericalContinuation(Week17Functions.ModifiedBetaFormHopf,[2.3,0,30],2,ParamNSteps=2,ParamStepSize=-0.1 )
    plt.plot(Space,Answer[:,0])
    plt.show()
    return       



if __name__ == "__main__":
    Main()