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
    """ Finds the finite difference solution to a diffusion type PDE 
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
    if not (Guess is None): #Checking for an approximate solution
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
    XSpace = XValues         # Cropping the range of x to match BCs
    if LeftBC[0][0]=="D":
        XSpace = XSpace[1:]
        U = U[1:]
    if RightBC[0][0]=="D":
        XSpace = XSpace[:-1]  
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

