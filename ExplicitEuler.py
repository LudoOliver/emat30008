"""Week20: Explicit Euler"""
from CommonModules import *
from EulerStep import EulerStep
import FiniteDifferences

import matplotlib.animation as animation
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
                                                    NPoints=NPoints+1,
                                                    Robin=Robin,
                                                    Neuman=Neuman,
                                                    Wrapped=1)
    
    XPoints = np.linspace(X0,XN,num=NPoints)
    TimePoints = np.linspace(0,TimeLimit,num=NTimeSteps)
    InitialU = InitalConditions(XPoints)
    #if InitialU[0]!=U0 or InitialU[-1]!=UN
    print(NTimeSteps)
    ULines = np.empty((NTimeSteps,NPoints))
    print(np.shape(ULines))
    ULines[0,:] = InitialU
    for i in range(1,NTimeSteps):
        ULines[i,:] = ULines[i-1,:]+FwdCoefficient*(Adx.dot(ULines[i-1,:])+Bdx)
        
    return XPoints,TimePoints,ULines

def SinICs(x,a,b):
    return np.sin(math.pi*(x-a)/(b-a))

def ExpectedSoln(x,t,a,b,d):
    #SolnLine = np.empty((np.size(x),np.size(t)))
    Xline = np.sin(math.pi*(x-a)/(b-a))
    #for i in range(0,len(t)):
    Tline = np.exp(-d*(math.pi**2)*t/((b-a)**2))
    SolnLine = np.outer(Xline,Tline).T
    return SolnLine
if __name__=="__main__":
    X0 = 0
    U0 = 0
    XN = 1
    UN = 0
    
    delta =1
    gamma =1
    
    DiffusionConstant = 1
    DlechtBCs = np.array([[X0,XN],[U0,UN]]) 
    NeumanBCs = np.array([[X0,XN],[U0,delta]]) #For du/dx|Xn = delta
    RobinBCS = np.array([[X0,XN],[U0,delta],[0,gamma]]) 
    
    SinICsAtBcs = lambda x: SinICs(x,a=X0,b=XN)
    X,T,Soln = ExplicitEuler(DlechtBCs
                            ,DiffusionConstant=DiffusionConstant
                            ,InitalConditions=SinICsAtBcs)
    TrueSoln = ExpectedSoln(X,T,a=X0,b=XN,d=DiffusionConstant)
    
    num_columns = 100
    num_rows = 40000
    def Animate(i):
        plt.plot(Soln[i,:])
        plt.plot(TrueSoln[i,:]+0.1)
    fig = plt.figure()
    ax = plt.axes(xlim=(0, 1), ylim=(0, 1))
    line1, = ax.plot([], [], label="Soln", lw=2)
    line2, = ax.plot([], [], label="True Soln", lw=2)
    #ax.xlim((0,1))
    #ax.ylim((0,1))
    def update(frame):
        line1.set_data(X, Soln[frame*200,:])
        line2.set_data(X, TrueSoln[frame*200,:])
        return line1, line2
    ani = animation.FuncAnimation(fig, update, frames=int(num_rows/200), blit=True)
    ax.legend()
    plt.show()