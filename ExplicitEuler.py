"""Week20: Explicit Euler"""
from CommonModules import *
from EulerStep import EulerStep
import FiniteDifferences
import ODESolver

import matplotlib.animation as animation

def MethodOfLines(BoundaryConditions,
                    Euler = 1,
                    Runge = 0,
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
    
    DeltaX = abs(XN-X0)/NPoints
    
    
    if NTimeSteps is None:
        DeltaT = (1/DiffusionConstant)*(DeltaX**2)*(1/2)
        NTimeSteps = math.ceil(TimeLimit/DeltaT)
    else:
        DeltaT = TimeLimit/NTimeSteps
    
    FwdCoefficient = DiffusionConstant*DeltaT/(DeltaX**2)
    #print(f"Delta t is {DeltaT}")
    if FwdCoefficient > 1/2:
        print(f"\n Error: Please modify your step values\n")
        return 0
    
    Adx, Bdx = FiniteDifferences.FiniteSolvePoisson(Bounds=BoundaryConditions,
                                                    NPoints=NPoints+1,
                                                    Robin=Robin,
                                                    Neuman=Neuman,
                                                    Wrapped=1)
    
    if Robin:
        0
    elif Neuman:
        0
    else:
        UN = Bdx[-1]
    XPoints = np.linspace(X0,XN,num=NPoints)
    TimePoints = np.linspace(0,TimeLimit,num=NTimeSteps)
    InitialU = InitalConditions(XPoints)
    #if InitialU[0]!=U0 or InitialU[-1]!=UN
    #print(NTimeSteps)
    ULines = np.empty((NTimeSteps,NPoints))
    if Euler :
        #print(np.shape(ULines))
        ULines[0,:] = InitialU
        for i in range(1,NTimeSteps):
            ULines[i,:] =  ULines[i-1,:]+(FwdCoefficient*(Adx.dot(ULines[i-1,:])+Bdx))+(DeltaX**2)*ReactionTerm(XPoints,ULines[i-1,:])
            ULines[i,0]=U0
            ULines[i,-1]=UN
    elif Runge :
        System = lambda u,t : (DiffusionConstant/(DeltaX**2))*(Adx.dot(u)+Bdx)+(1/FwdCoefficient)*ReactionTerm(XPoints,u)
        ULines,TimePoints = ODESolver.Solve_to(System,
                                    x0=InitialU,
                                    tspan=(0,TimeLimit),
                                    DeltaTMax=DeltaT
                                    )
        ULines = ULines.T 
    else:
        print(f"Please Specify Solver")
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

def Main1():
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
    X,T,Soln = MethodOfLines(DlechtBCs
                            ,DiffusionConstant=DiffusionConstant
                            ,InitalConditions=SinICsAtBcs)
    TrueSoln = ExpectedSoln(X,T,a=X0,b=XN,d=DiffusionConstant)
    
    Xrk,Trk,Solnrk = MethodOfLines(DlechtBCs
                            ,DiffusionConstant=DiffusionConstant
                            ,InitalConditions=SinICsAtBcs,Euler=0,Runge=1)

    num_columns = 100
    num_rows = 40000

    fig = plt.figure()
    ax = plt.axes(xlim=(0, 1), ylim=(0, 1))
    line1, = ax.plot([], [], label="Soln", lw=2)
    line2, = ax.plot([], [], label="True Soln", lw=2)
    line3, = ax.plot([],[],label="Soln with runge-kutta",lw=2)
    #ax.xlim((0,1))
    #ax.ylim((0,1))
    def update(frame):
        line1.set_data(X, Soln[frame*200,:])
        line2.set_data(X, TrueSoln[frame*200,:])
        line3.set_data(X, Solnrk[frame*200,:])
        return line1, line2, line3
    ani = animation.FuncAnimation(fig, update, frames=int(num_rows/200), blit=True)
    ax.legend()
    plt.show()

def BratuICs(x):
    return 0*x

def BratuSource(x,u,mu=0.1):
    y = np.exp(u*mu)#,dtype=np.longdouble)
    #if np.isnan(y).any():
        #print(f"Failed at {u.max()}")
        #return 
    return y
def Main2(): #Trying
    X0 = 0
    U0 = 0
    XN = 1
    UN = 0
    
    DiffusionConstant = 1
    DlechtBCs = np.array([[X0,XN],[U0,UN]]) 
    
    X,T,Soln = MethodOfLines(DlechtBCs,
                            ReactionTerm=BratuSource,
                            InitalConditions=BratuICs)
    #print(Soln)
    plt.plot(X,Soln[-1,:])
    #plt.plot(X,Soln[40,:])
    plt.title("Euler")
    plt.show()
def Main3():
    X,T,Soln = MethodOfLines(DlechtBCs,
                            ReactionTerm=BratuSource,
                            InitalConditions=BratuICs,
                            Euler=0,
                            Runge=1)
    #print(Soln.max())
    plt.plot(X,Soln[-1,:])
    plt.title("Runge")
    plt.show()
if __name__=="__main__":
    #Main1()
    X0 = 0
    U0 = 0
    XN = 1
    UN = 0
    
    DiffusionConstant = 1
    DlechtBCs = np.array([[X0,XN],[U0,UN]]) 
    Main3()
    Main2()