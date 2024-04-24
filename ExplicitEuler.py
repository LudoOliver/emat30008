"""Week20: Explicit Euler"""
import scipy.optimize
from CommonModules import *
from EulerStep import EulerStep
import FiniteDifferences
import ODESolver
from AMatrixBVector import MakeAMatrixBVector
import matplotlib.animation as animation

def MethodOfLines(  LeftBC,LeftBCLocation,
                    RightBC, RightBCLocation,
                    Euler = 1,
                    Runge = 0,
                    DiffusionConstant =1,
                    TimeLimit = 1,
                    ReactionTerm = FiniteDifferences.NoSourceTerm,
                    NTimeSteps = None,
                    InitalConditions = None,
                    Implicit = 0,
                    NPoints = 100,
                    ):
    """ Uses numerical methods to evaluate a PDE over time
    Parameters
        [Left/Right]Bc (list): list of form (Type,(Parameters))
                Type: Either "D","N" or"R" for Dirlecht,Neuman or Robin respectively
                Parameters: scalar for "D" or "N", tuple (delta,gamma) for "R"
        [Left/Right]Location (scalar):
            Position of the left/right boundary in x domain
        DiffusionConstant (scalar): D - default is 1
        Euler(boolean,optional) : Euler method on/off, on by default
        Implicit(boolean,optional) : the choice to use implicit methods, off by default
        Runge(boolean,optional) : the choice to use RK4 instead - Euler must be turned off to use
        Reaction (f(x,u),optional) : Reaction term - off by default
        Guess (array,optional) : Solution estimate 
        NPoints(int,optional): Number of points in matrix, default is 100
    Returns:
        X(array): x values of each grid point in the domain
        U(array): U(x) for  x in X
    """
    DeltaX = (RightBCLocation-LeftBCLocation)/NPoints
    
    if NTimeSteps is None:
        DeltaT = (1/DiffusionConstant)*(DeltaX**2)*(1/2)
        NTimeSteps = math.ceil(TimeLimit/DeltaT)
    else:
        DeltaT = TimeLimit/NTimeSteps
    
    FwdCoefficient = DiffusionConstant*DeltaT/(DeltaX**2) 
    #print(f"Delta t is {DeltaT}")
    if FwdCoefficient > 1/2 and not Implicit:
        print(f"\n Error: Please modify your step values\n")
        return 0
    

    AMatrix, BVector = MakeAMatrixBVector(NPoints=NPoints,DeltaX=DeltaX,
                                            Left=LeftBC,Right=RightBC) 

    
    XPoints = np.linspace(LeftBCLocation,RightBCLocation,num=NPoints)
    TimePoints = np.linspace(0,TimeLimit,num=NTimeSteps)
    InitialU = InitalConditions(XPoints)

    XSpace = XPoints
    if LeftBC[0][0]=="D":   ##Ensuring correct Grid dimensions for the boundary conditions
        XSpace = XSpace[1:]
        InitialU = InitialU[1:]
    if RightBC[0][0]=="D":
        XSpace = XSpace[:-1]  
        InitialU = InitialU[:-1]

    ULines = np.empty((NTimeSteps,len(XSpace))) #Initialising the array to store solutions
    
    
    if Euler :
        ULines[0,:] = InitialU
        for i in range(1,NTimeSteps):
            if Implicit:
                FImplicit = lambda Soln : ((Soln-ULines[i-1,:])/DeltaT)-(DiffusionConstant/(DeltaX**2))*(AMatrix.dot(Soln)+BVector)-DeltaT*ReactionTerm(XSpace,Soln)
                ImplicitSoln = scipy.optimize.root(FImplicit,ULines[i-1,:])
                if ImplicitSoln.success:  #Implicit Euler Method, with additional error handling options 
                    ULines[i,:] = ImplicitSoln.x
                else:
                    if LeftBC[0][0]=="D":
                                ULines = np.hstack((LeftBC[1]*np.ones(NTimeSteps)[:,None],ULines))
                    if RightBC[0][0]=="D":
                                ULines = np.hstack((ULines,RightBC[1]*np.ones(NTimeSteps)[:,None]))
                    print(f" \nFailed:\n   {ImplicitSoln.message} \n   Returned Values for {i} steps \n")
                    return XPoints[:i],TimePoints[:i],ULines[:i,:]
            else:  #Explicit Euler
                ULines[i,:] =  ULines[i-1,:]+(FwdCoefficient*(AMatrix.dot(ULines[i-1,:])+BVector))+(DeltaX**2)*ReactionTerm(XSpace,ULines[i-1,:]) #Might want delta x**2

        if LeftBC[0][0]=="D":   #Forcing Dirilecht boundaries
            ULines = np.hstack((LeftBC[1]*np.ones(NTimeSteps)[:,None],ULines))
        if RightBC[0][0]=="D":
            ULines = np.hstack((ULines,RightBC[1]*np.ones(NTimeSteps)[:,None]))
    elif Runge :  #Implementing the Runge Kutta solution
        System = lambda u,t : (DiffusionConstant/(DeltaX**2))*(AMatrix.dot(u)+BVector)+ReactionTerm(XPoints,u)
        ULines,TimePoints = ODESolver.Solve_to(System,
                                    x0=InitialU,
                                    tspan=(0,TimeLimit),
                                    DeltaTMax=DeltaT
                                    )
        ULines = ULines.T 
    else:
        print(f"Please Specify Solver")  #Ensuring a solver was specified
    return XPoints,TimePoints,ULines






