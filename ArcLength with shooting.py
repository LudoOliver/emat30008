from CommonModules import *
import Week16General
import Week17Functions
import ODESolver

def ShootingSolveWrapper(Func,t,Param,SolnEstimate):
    #Not intended for front end use
    WrappedFunc = lambda u,t :Func(u, t,Param)
    #print(SolnEstimate[:-1])
    #print(SolnEstimate[-1])
    print(Param)
    #print(SolnEstimate)
    x,period = Week16General.Shooting(WrappedFunc, SolnEstimate[:-1], SolnEstimate[-1])#,StepSize=0.0001)
    return np.hstack([x,period])


def ShootingArcLengthCont(Func,X0,ParamBounds,ContinuationMaxSteps,
                        ParamStepSize=0.1,Solver=ODESolver.RungeKutta4,
                        SolverStepSize=0.1
                        ):
    
    """ Func(function) : function of (u,x)
        X0 (array): approximate  initial conditions of form (x1,..xn-1,T)
        ParamBounds (tuple):the region of parameter space to investigate
        ContinuationMaxSteps (int): the upper bound on the number of steps taken
        ParamStepSize (float,optional): initial step in parameter space, default is 0.1
        Solver (func,optional) : ode solver to use, default is RK4
        SolverStepSize (float,optional) : step size to solve ODEs with, default is 0.1
    Returns:
        SolnSpace (array): 2d array of solutions
        ParamSpace (array): 1d array of the corresponding parameter values
    """
    
    P0,PN = ParamBounds[0],ParamBounds[1]
    Direction = np.sign(PN-P0)
    #ParamArray = np.linspace(P0,PN,num=MaxSteps)
    ParamSpace = np.empty(ContinuationMaxSteps)+np.nan
    SolnSpace = np.empty([ContinuationMaxSteps,np.size(X0)])+np.nan
    
    ParamSpace[0],ParamSpace[1] = ParamBounds[0],ParamBounds[0]+Direction*ParamStepSize
    SolnSpace[0,:] = ShootingSolveWrapper(Func,1,ParamSpace[0],X0)
    SolnSpace[1,:] = ShootingSolveWrapper(Func,1,ParamSpace[1],SolnSpace[0,:])

    def ArcShootRootFind(SolnAndParam):
        
        SolnToInvestigate,ParamToInvestigate = SolnAndParam[:-1],SolnAndParam[-1] 
        #Unpacking u,P array - paired together for ease of use with scipy solve
        FuncAtParam = lambda x,t : Func(x,t,ParamToInvestigate)
        #SingleShotArgs = (SolverStepSize,Solver,FuncAtParam)
        #Fixing function to current param value
        ShootingLHS = Week16General.SingleShot(SolnToInvestigate,SolverStepSize,Solver,FuncAtParam)
        #Shooting root finding problem
        ArcCondition = (np.dot(ParamSecant,(ParamToInvestigate-ParamEstimate))+
                                np.dot(SolnSecant,(SolnToInvestigate-SolnEstimate)))
        #Arc length condition, using previous solutions from the local scope
        ArcLengthRHS = np.append(ShootingLHS,ArcCondition)
        
        return ArcLengthRHS 
    
    for i in range(2,ContinuationMaxSteps):
        if min(ParamBounds) < ParamSpace[i-1] < max(ParamBounds):
            ParamSecant = ParamSpace[i-1]-ParamSpace[i-2]
            SolnSecant = SolnSpace[i-1,:]-SolnSpace[i-2,:]
            SolnEstimate = SolnSpace[i-1,:]+SolnSecant
            ParamEstimate = ParamSpace[i-1]+ParamSecant 
            #Using the scope of nested function definitions to reduce difficulty passing arguments

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
            #If an error is encountered, all results from up to that point are returned
            return SolnSpace,ParamSpace 
            
    SolnSpace = np.squeeze([i for i in SolnSpace if not np.isnan(i).all()]) #Ensuring correct output dimensions
    ParamSpace = [i for i in ParamSpace if not np.isnan(i)]    
    return SolnSpace,ParamSpace

