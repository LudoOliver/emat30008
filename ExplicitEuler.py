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
    
    # Adx, Bdx = FiniteDifferences.FiniteSolvePoisson(Bounds=BoundaryConditions,
    #                                                NPoints=NPoints+1,
    # #                                               Robin=Robin,
    #  #                                              Neuman=Neuman,
    #                                             Wrapped=1)
    AMatrix, BVector = MakeAMatrixBVector(NPoints=NPoints,DeltaX=DeltaX,
                                            Left=LeftBC,Right=RightBC)
    #NPoints = len(BVector) # Leaving changing Grid size to old function
    
    XPoints = np.linspace(LeftBCLocation,RightBCLocation,num=NPoints)
    TimePoints = np.linspace(0,TimeLimit,num=NTimeSteps)
    InitialU = InitalConditions(XPoints)
    #XValues = np.linspace(LeftBCLocation,RightBCLocation, num=NPoints)
    XSpace = XPoints
    if LeftBC[0][0]=="D":
        XSpace = XSpace[1:]
        InitialU = InitialU[1:]
    if RightBC[0][0]=="D":
        XSpace = XSpace[:-1]  #Cant index end with 0 nicely otherwise could do inline
        InitialU = InitialU[:-1]
    #if InitialU[0]!=U0 or InitialU[-1]!=UN
    #print(NTimeSteps)
    ULines = np.empty((NTimeSteps,len(XSpace)))
    
    
    if Euler :
        #print(np.shape(ULines))
        ULines[0,:] = InitialU
        for i in range(1,NTimeSteps):
            if Implicit:
                FImplicit = lambda Soln : ((Soln-ULines[i-1,:])/DeltaT)-(DiffusionConstant/(DeltaX**2))*(AMatrix.dot(Soln)+BVector)-DeltaT*ReactionTerm(XSpace,Soln)
                ImplicitSoln = scipy.optimize.root(FImplicit,ULines[i-1,:])
                if ImplicitSoln.success:
                    ULines[i,:] = ImplicitSoln.x
                else:
                    if LeftBC[0][0]=="D":
                                ULines = np.hstack((LeftBC[1]*np.ones(NTimeSteps)[:,None],ULines))
                    if RightBC[0][0]=="D":
                                ULines = np.hstack((ULines,RightBC[1]*np.ones(NTimeSteps)[:,None]))
                    print(f" \nFailed:\n   {ImplicitSoln.message} \n   Returned Values for {i} steps \n")
                    return XPoints[:i],TimePoints[:i],ULines[:i,:]
        #print(f")
            else:
                ULines[i,:] =  ULines[i-1,:]+(FwdCoefficient*(AMatrix.dot(ULines[i-1,:])+BVector))+(DeltaX**2)*ReactionTerm(XSpace,ULines[i-1,:]) #Might want delta x**2
            #ULines[i,0]=U0
            #ULines[i,-1]=UN
        if LeftBC[0][0]=="D":
            ULines = np.hstack((LeftBC[1]*np.ones(NTimeSteps)[:,None],ULines))
        if RightBC[0][0]=="D":
            ULines = np.hstack((ULines,RightBC[1]*np.ones(NTimeSteps)[:,None]))
    elif Runge :
        System = lambda u,t : (DiffusionConstant/(DeltaX**2))*(AMatrix.dot(u)+BVector)+ReactionTerm(XPoints,u)
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
    
    NewBcs = ("D",0)
    
    
    SinICsAtBcs = lambda x: SinICs(x,a=X0,b=XN)
    X,T,Soln = MethodOfLines(LeftBC=NewBcs,LeftBCLocation=0,
                            RightBC=NewBcs,RightBCLocation=1,
                            DiffusionConstant=DiffusionConstant,
                            Runge=1,
                            InitalConditions=SinICsAtBcs)
    
    TrueSoln = ExpectedSoln(X,T,a=X0,b=XN,d=DiffusionConstant)
    
    ImpX,ImpT,ImpSoln = MethodOfLines(LeftBC=NewBcs,LeftBCLocation=0,
                            RightBC=NewBcs,RightBCLocation=1,
                            DiffusionConstant=DiffusionConstant,
                            Implicit=1,
                            InitalConditions=SinICsAtBcs)

    num_columns = 100
    num_rows = 40000
    plt.figure()
    plt.plot(ImpX,ImpSoln[-1,:])
    
    fig = plt.figure()
    ax = plt.axes(xlim=(0, 1), ylim=(0, 1))
    line1, = ax.plot([], [], label="Soln", lw=2)
    line2, = ax.plot([], [], label="True Soln", lw=2)
    #line3, = ax.plot([],[],label="Implicit Soln",lw=2)
    #ax.xlim((0,1))
    #ax.ylim((0,1))
    def update(frame):
        line1.set_data(X, Soln[frame*200,:])
        line2.set_data(X, TrueSoln[frame*200,:])
        #line3.set_data(ImpX, ImpSoln[frame*200,:])
        return line1, line2#, line3
    ani = animation.FuncAnimation(fig, update, frames=int(num_rows/200)-1, blit=True)
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
    NewBcs = ("D",0)
    X,T,Soln = MethodOfLines(LeftBC=NewBcs,LeftBCLocation=0,
                            RightBC=NewBcs,RightBCLocation=1,
                            ReactionTerm=BratuSource,
                            InitalConditions=BratuICs)
    #print(Soln)
    plt.plot(X,Soln[-10,:])
    plt.plot(X,Soln[-1,:])
    plt.xlim(0,1)
    #plt.plot(X,Soln[40,:])
    plt.title("Euler")
    plt.show()
def Main3():
    NewBcs = ("D",0)
    X0 = 0
    U0 = 0
    XN = 1
    UN = 0
    SinICsAtBcs = lambda x: SinICs(x,a=X0,b=XN)
    DiffusionConstant = 1
    ImpX,ImpT,ImpSoln = MethodOfLines(LeftBC=NewBcs,LeftBCLocation=0,
                            RightBC=NewBcs,RightBCLocation=1,
                            DiffusionConstant=DiffusionConstant,
                            NTimeSteps=100,
                            Implicit=1,
                            InitalConditions=SinICsAtBcs)

    print(len(ImpX))
    plt.figure()
    for i in range(0,len(ImpX)):
        if True:#i%(len(ImpX)/30)==0:
            plt.plot(ImpX,ImpSoln[i,:],label=str(i))
            plt.xlim(0,1)
    plt.legend()
    plt.show()
if __name__=="__main__":
    Main1()
    #X0 = 0
    #U0 = 0
    #XN = 1
    #UN = 0
    
    #DiffusionConstant = 1
    #DlechtBCs = np.array([[X0,XN],[U0,UN]]) 
    #Main3()
    #Main2()