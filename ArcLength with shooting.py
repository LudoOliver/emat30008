from CommonModules import *
import Week16General
import Week17Functions
import ODESolver

def ShootingSolveWrapper(Func,t,Param,SolnEstimate):
    WrappedFunc = lambda u,t :Func(u, t,Param)
    #print(SolnEstimate[:-1])
    #print(SolnEstimate[-1])
    print(Param)
    #print(SolnEstimate)
    x,period = Week16General.Shooting(WrappedFunc, SolnEstimate[:-1], SolnEstimate[-1])#,StepSize=0.0001)
    return np.hstack([x,period])


def ShootingArcLengthCont(Func,X0,ParamBounds,ContinuationMaxSteps,
                        ParamStepSize=0.1,WithShooting=1,Solver=ODESolver.RungeKutta4,
                        SolverStepSize=0.1
                        ):
    
    """_summary_

    Returns:
        _type_: _description_
    """
    
    P0,PN = ParamBounds[0],ParamBounds[1]
    Direction = np.sign(PN-P0)
    #ParamArray = np.linspace(P0,PN,num=MaxSteps)
    ParamSpace = np.empty(ContinuationMaxSteps)+np.nan
    SolnSpace = np.empty([ContinuationMaxSteps,np.size(X0)])+np.nan
    
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
        
        SolnToInvestigate,ParamToInvestigate = SolnAndParam[:-1],SolnAndParam[-1]
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
            #print([SolnEstimate,ParamEstimate])
            #print(np.append(SolnEstimate,ParamEstimate))
            SolnAndParam = np.append(SolnEstimate,ParamEstimate)
            Result = scipy.optimize.root(ArcShootRootFind,SolnAndParam).x
            print(Result)
            #NewSoln,NewParam = Result[0],Result[1]
            SolnSpace[i,:] = Result[:-1]
            ParamSpace[i] = Result[-1]
            
        else:
            SolnSpace = np.squeeze([i for i in SolnSpace if not np.isnan(i).all()])
            ParamSpace = [i for i in ParamSpace if not np.isnan(i)]
            print("Bounds Exceeded")
            return SolnSpace,ParamSpace
    SolnSpace = np.squeeze([i for i in SolnSpace if not np.isnan(i).all()])
    ParamSpace = [i for i in ParamSpace if not np.isnan(i)]    
    return SolnSpace,ParamSpace

def Main(): 
    SolnSpace,ParamSpace = ShootingArcLengthCont(Week17Functions.BetaFormHopf,[1.1,0.1,30],[0.4,2],15,SolverStepSize=0.001)
    #ShootingNumericalContinuation(Week17Functions.ModifiedBetaFormHopf,[2.3,0,30],2,ParamNSteps=2,ParamStepSize=-0.1 )
    plt.plot(ParamSpace,SolnSpace[:,0])
    plt.show()
    return    
def TestWrapper():
    Ans1,Time1 = Week16General.Shooting(Week17Functions.BetaFormHopf,[1.1,0],30)#,StepSize=0.0000001)
    print(Ans1)
    print(Time1)
    Ans = ShootingSolveWrapper(Week17Functions.BetaFormHopf,1,0,[1.1,0.1,6])
    print(Ans)
    print(Ans1)
if __name__ == "__main__":
    #TestWrapper()
    Main()