# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 20:26:26 2024

@author: Admin
"""

from CommonModules import *
import ODESolver
from AMatrixBVector import MakeAMatrixBVector
def FirstOrderFiniteDiff(Func,StepSize,LocationPoint):
    dx = (Func(LocationPoint+StepSize)-Func(LocationPoint-StepSize))
    dt = (2*StepSize)
    return dx/dt
def SecondOrderFiniteDiff(Func,StepSize,LocationPoint):
    ddx = Func(LocationPoint+StepSize)-2*Func(LocationPoint)+Func(LocationPoint-StepSize)
    dtt = StepSize**2
    return ddx/dtt

def NoSourceTerm(x,U):
    return 0

def SimpleSourceTerm(x,U):
    return 1 # Probably could just be 

def BratuTerm(x,U,mu=2):
    return np.exp(mu*U)
#This is no longer in use
def FiniteSolvePoisson(Bounds, 
                        DiffusionConstant=1,Reaction=NoSourceTerm,
                        NPoints=100,Guess=None,
                        Wrapped=0,
                        Robin=0,Neuman=0):
    """_summary_

    Args:
        Bounds (2d array): boundary conditions of diffusion equation
        DiffusionConstant (int, optional): _description_. Defaults to 1.
        Reaction (f(x,U), optional): reaction term
        NPoints (int, optional): grid size Defaults to 100.
        Guess (Array, optional): estimated array of solutions
        Wrapped (binary, optional): output only A and B arrays - only for internal use
        Robin (binary, optional): Use Robin BCs
        Neuman (binary, optional): Use Neuman BCs

    Returns:
        X(array) : length n vector of x co-ordinates
        U(array) : length n vector of u(x)
        
    if Wrapped:
        Returns: 
            A^d matrix
            B^d vector
    """
    
    X0 = Bounds[0,0]
    U0 = Bounds[1,0]
    XN = Bounds[0,1]
    if Guess is None: # Slightly inconsistent but the best solution as numpy arrays dont have a set truth value   

        DeltaX = (XN-X0)/NPoints
        U = np.zeros(NPoints)
        
    else:
        #NPoints = math.ceil((abs(XN-X0))/MaxStepSize)
        #DeltaX = (XN-X0)/NPoints
        #U = np.zeros(NPoints+1)
        NPoints = len(Guess)
        DeltaX = (XN-X0)/NPoints
        U = Guess
        
    #TrueStepSize = (XN-X0)/NPoints
    
    if Neuman and Robin:
        print("Error: Please specify only one type of boundary conditions")
        return 0
    elif Robin:
        delta = Bounds[1,1]
        UN = 2*delta*DeltaX
        gamma = Bounds[2,1]
        EndVector = (2,-2*(1+gamma*DeltaX))
    elif Neuman:
        delta = Bounds[1,1]
        UN = 2*delta*DeltaX
        EndVector = (2,-2)
    else:
        NPoints=NPoints-1
        U= U[:-1]
        UN = Bounds[1,1]
        EndVector = (1,-2)    
    #print(U)
    XSpace = np.linspace(X0,XN,num=NPoints)
    #print("Size",np.size(XSpace))
    #return
    #print(XSpace[2]-XSpace[1]+TrueStepSize)
    CoeffMatrix = np.zeros([NPoints,NPoints])
    CoeffMatrix[0,(0,1)] = -2,1
    CoeffVector = np.array([-1,0,1])
    for i in range(1,len(CoeffMatrix)-1):
        CoeffMatrix[i,i+CoeffVector] = (1,-2,1)
    CoeffMatrix[-1,(-2,-1)] = EndVector
    #print(f"Grid length : {np.sqrt(np.size(CoeffMatrix))}")
    
    #Maybe u=zeros
    U[0] = U0
    U[-1] = UN
    
    if Wrapped:
        return CoeffMatrix,U
        
    
    def EqnToSolve(SolnArray):
        return DiffusionConstant*(CoeffMatrix.dot(SolnArray)+U)+Reaction(XSpace,SolnArray)*(DeltaX**2) #Working
        #return DiffusionConstant*(1/DeltaX**2)*(CoeffMatrix.dot(SolnArray)+U)+Reaction(XSpace,SolnArray)
    Ans = scipy.optimize.root(EqnToSolve,U)
    
    #Ans = Ans
    #print(Ans[0],Ans[1],Ans[0]-U[0],Ans[1]-Ans[0])
    #print(Ans[-1]+(Ans[1]-Ans[0]))
    if Ans.success:
        return XSpace, Ans.x
    else:
        print(f" \nFailed:\n   {Ans.message} \n   Try improving BCs\n")
        #print(f")
        return 0,0

def ExpectedSoln(BoundaryCond,Xarray,DiffusionConst=1):
    X0 = BoundaryCond[0,0]
    XN = BoundaryCond[0,1]
    U0 = BoundaryCond[1,0]
    UN = BoundaryCond[1,1]
    Soln = np.add((-1/(2*DiffusionConst))*np.multiply((Xarray-X0),(Xarray-XN)),(UN-U0)/(XN-X0)*(Xarray-X0)+U0)
    return Soln
#%%
def FiniteDifferences(  LeftBC,LeftBCLocation,
                        RightBC, RightBCLocation,
                        DiffusionConstant=1,Reaction=NoSourceTerm,

                        NPoints=100,Guess=None):
    """
    Parameters
        [Left/Right]Bc (list): list of form (Type,(Parameters))
                Type: Either "D","N" or"R" for Dirlecht,Neuman or Robin respectively
                Paramters: scalar for "D" or "N", tuple (delta,gamma) for "R"
        [Left/Right]Location (scalar):
            Postition of left/right boundary in x domain
        DiffusionConstant (scalar): D - default is 1
        Reaction (f(x,u),optional) : Reaction term - off by default
        Guess (array,optional) : Solution estimate 
        NPoints(int,optional): Number of points in matrix, default is 100
    Returns:
        X(array): x values of each grid point in the domain
        U(array): U(x) for  x in X
    """
    if not (Guess is None): 
        NPoints = len(Guess)
        U = Guess
        UFromGuess = 1
    else:
        UFromGuess = 0
        U = np.zeros(NPoints)
    
    DeltaX = (RightBCLocation-LeftBCLocation)/NPoints
    AMatrix , BVector = MakeAMatrixBVector(NPoints=NPoints,DeltaX=DeltaX,
                                            Left=LeftBC,Right=RightBC,
                                            FromGuess=UFromGuess)
    
    XValues = np.linspace(LeftBCLocation,RightBCLocation, num=NPoints)
    XSpace = XValues
    if LeftBC[0][0]=="D":
        XSpace = XSpace[1:]
        U = U[1:]
    if RightBC[0][0]=="D":
        XSpace = XSpace[:-1]  #Cant index end with 0 nicely otherwise could do inline
        U = U[:-1]


    def MatrixSystem(Soln):
            return DiffusionConstant*(AMatrix.dot(Soln)+BVector)+Reaction(XSpace,Soln)*(DeltaX**2)
    
    SolnU = scipy.optimize.root(MatrixSystem,U)
    if not SolnU.success:
        #return XSpace, SolnU.x
    # else:
        print(f" \nFailed:\n   {SolnU.message} \n   Check Bcs \n")
        #print(f")
        return 0,0
    
    U = SolnU.x
    if LeftBC[0][0]=="D":
        U = np.r_[LeftBC[1],U]
    if RightBC[0][0]=="D":
        U = np.r_[U,RightBC[1]]
    
    return XValues, U

#%%
if __name__ == "__main__":
    
    
    NewBounds = ("D",(0))
    NewX,NewU = FiniteDifferences(LeftBC=NewBounds,LeftBCLocation=0,
                                    RightBC=NewBounds,RightBCLocation=1,
                                    Reaction=BratuTerm)
    


    plt.figure()
    plt.plot(NewX,NewU,label="New soln")

    plt.legend()
    plt.show()