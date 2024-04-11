"""Week20: Explicit Euler"""
from CommonModules import *
from EulerStep import EulerStep
import FiniteDifferences
def ExplicitEuler(BoundaryConditions,
                    DiffusionConstant =1,
                    TimeLimit = 1,
                    ReactionTerm = FiniteDifferences.NoSourceTerm,
                    NTimeSteps = None,
                    InitalConditions = None,
                    NPoints = 100,
                    Robin=0,Neuman=0
                    ):
    X0 = BoundaryConditions[0,0]
    U0 = BoundaryConditions[1,0]
    XN = BoundaryConditions[0,1]
    DeltaX = (X0-XN)/NPoints
    if NTimeSteps is None:
        DeltaT = (1/DiffusionConstant)*(DeltaX**2)*(1/4)
        NTimeSteps = math.ceil(TimeLimit/DeltaT)
    else:
        DeltaT = TimeLimit/NTimeSteps
    
    FwdCoefficient = DiffusionConstant*DeltaT/(DeltaX**2)
    
    if FwdCoefficient > 1/2:
        print(f"\n Error: Please modify your step values\n")
        return 0
    
    Adx, Bdx = FiniteDifferences.FiniteSolvePoisson(Bounds=BoundaryConditions,
                                                    NPoints=NPoints,
                                                    Robin=Robin,
                                                    Neuman=Neuman,
                                                    Wrapped=1)
    
    XPoints = np.linspace(X0,XN,num=NPoints)
    InitialU = InitalConditions(XPoints)
    #if InitialU[0]!=U0 or InitialU[-1]!=UN
    ULines = np.empty((NPoints,NTimeSteps))
    ULines[0,:] = InitialU
    for i in range(1,NTimeSteps+1):
        ULines[i,:] = ULines[i-1,:]+FwdCoefficient*(Adx.dot(ULines[i-1,:])+Bdx)
        
    return XPoints, InitialU

def SinICs(x,a,b):
    return np.sin(math.pi*(x-a)/(b-a))