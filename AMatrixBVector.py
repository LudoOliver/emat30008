import numpy as np 

def MakeAMatrixBVector(NPoints,DeltaX,Left,Right,FromGuess=0):
    """Create sthe Axx matrix and Bxx vector needed for finite difference methods

    Args:
        NPoints (int): number of points used to define the gird
        DeltaX (float): step size in x
        Left (tuple): left boundary condition
        Right (tuple): right boundary condition_
        FromGuess (int, optional): _description_. Defaults to 0.

    Returns:
        A (2d-array): the Axx matrix 
        B (1d-array): the Bxx vector
    """
    if Left[0][0] == "R":  #Double index allows spelling errors, using properties of length(1) strs
        D1 = Left[1][0]
        Gamma1 = Left[1][1]
        ALeftCorner = (-2*(1-Gamma1*DeltaX),2)
        B0 = -2*D1*DeltaX
    elif Left[0][0] == "N":
        D1 = Left[1]
        ALeftCorner = (-2,2)
        B0 = -2*D1*DeltaX
    elif Left[0][0] == "D":
        NPoints = NPoints-1
        ALeftCorner = (-2,1)
        B0 = Left[1]
    else:
        print(f"\n Error: please check left BCs \n")
    
    if Right[0][0] == "R":
        D2 = Right[1][0]
        Gamma2 = Right[1][1]
        ARightCorner = (2,-2*(1-Gamma2*DeltaX))
        BN = 2*D2*DeltaX
    elif Right[0][0] == "N":
        D2 = Right[1]
        ARightCorner = (2,-2)
        BN = 2*D2*DeltaX
    elif Right[0][0] == "D":
        NPoints = NPoints-1
        ARightCorner = (1,-2)
        BN = Right[1]
        
    else:
        print(f"\n Error: please check right BCs \n")
    AMatrix = np.zeros([NPoints,NPoints])
    AMatrix[0,(0,1)] = ALeftCorner
    for i in range(1,NPoints-1):
        AMatrix[i,(i-1,i,i+1)]= (1,-2,1)
    AMatrix[-1,(-2,-1)] = ARightCorner
    
    BVector = np.zeros(NPoints)
    BVector[0] = B0
    BVector[-1] = BN
    
    return AMatrix, BVector

def Main():
    X0 = 0
    U0 = 0
    XN = 1
    UN = 0
    
    delta =1
    gamma =1
    
    DlechtBCs = np.array([[X0,XN],[U0,UN]]) 
    NeumanBCs = np.array([[X0,XN],[U0,delta]]) #For du/dx|Xn = delta
    RobinBCS = np.array([[X0,XN],[U0,delta],[0,gamma]]) 
    
    OldA,OldB = FiniteDifferences.FiniteSolvePoisson(Bounds=DlechtBCs,Wrapped=1)
    
    NewBCs = ("D",(0))
    NewA,NewB = MakeAMatrixBVector(99,0.01,Left=NewBCs,Right=NewBCs)
    print("New properties")
    print(f" A size = {np.size(NewA)}")
    print(f" B size = {np.size(NewB)}")
    print(f"\n Old Properties")
    print(f"A size = {np.size(OldA)}")
    print(f"B size = {np.size(OldB)}")
    print(f"\n Difference")
    print(f"A {np.sum(np.abs(OldA-NewA))}")
    print(f"B {np.sum(np.abs(OldB-NewB))} \n")
    print(OldB)
    print(NewB)
if __name__=="__main__":
    Main()