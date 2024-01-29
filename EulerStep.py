def EulerStep(Ftx,Xn,Tn,StepSize):
    return Xn+StepSize*Ftx(Xn,Tn),Tn+StepSize