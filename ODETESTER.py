"""Investigating ODEs"""
from CommonModules import *
from Week17Functions import *
import cProfile
import ODESolver
from CommonModules import *
import ExplicitEuler
LeftBC =("N",0)
RightBC = LeftBC
DiffusionConstant = 0.01
QuestionFiveICs = lambda x:0*x
def QuestionFiveSourceTerm(x,u):
    coef = np.square(1-u)
    exponent = np.exp(x*-1)
    return np.multiply(coef,exponent)

def Explicit():
    expX,expT,expU = ExplicitEuler.MethodOfLines(LeftBC=LeftBC,LeftBCLocation=0,
                                    RightBC=RightBC,RightBCLocation=6,
                                    NPoints=101,
                                    NTimeSteps=28000,
                                    TimeLimit=100,
                                    InitalConditions=QuestionFiveICs,
                                    ReactionTerm = QuestionFiveSourceTerm,
                                    DiffusionConstant=DiffusionConstant
                                    )
    return expX,expT,expU

def Runge():
    rungeX,rungeT,rungeU = ExplicitEuler.MethodOfLines(LeftBC=LeftBC,LeftBCLocation=0,
                                    RightBC=RightBC,RightBCLocation=6,
                                    NPoints=101,
                                    #NTimeSteps=1000, 
                                    TimeLimit=100,
                                    InitalConditions=QuestionFiveICs,
                                    ReactionTerm = QuestionFiveSourceTerm,
                                    DiffusionConstant=DiffusionConstant,
                                    Runge=1,
                                    Euler=0
                                    )
def Implicit():
    impX,impT,impU = ExplicitEuler.MethodOfLines(LeftBC=LeftBC,LeftBCLocation=0,
                                    RightBC=RightBC,RightBCLocation=6,
                                    NPoints=101,
                                    NTimeSteps=99,
                                    TimeLimit=100,
                                    InitalConditions=QuestionFiveICs,
                                    ReactionTerm = QuestionFiveSourceTerm,
                                    Implicit=1,
                                    DiffusionConstant=DiffusionConstant
                                    )
    return impX,impT,impU


# %%
print("imp")
cProfile.run("ExplicitEuler.MethodOfLines(LeftBC=LeftBC,LeftBCLocation=0,RightBC=RightBC,RightBCLocation=6,NPoints=101,NTimeSteps=99,TimeLimit=100,InitalConditions=QuestionFiveICs,ReactionTerm = QuestionFiveSourceTerm,Implicit=1,DiffusionConstant=DiffusionConstant)",sort="cumtime")

#%%
print("runge")
cProfile.run("ExplicitEuler.MethodOfLines(LeftBC=LeftBC,LeftBCLocation=0,RightBC=RightBC,RightBCLocation=6,NPoints=101,NTimeSteps=99,TimeLimit=100,InitalConditions=QuestionFiveICs,ReactionTerm = QuestionFiveSourceTerm,Implicit=1,DiffusionConstant=DiffusionConstant)",sort="cumtime")
                                    

#%%
print("exp")
cProfile.run("ExplicitEuler.MethodOfLines(LeftBC=LeftBC,LeftBCLocation=0,RightBC=RightBC,RightBCLocation=6,NPoints=101,NTimeSteps=28000,TimeLimit=100,InitalConditions=QuestionFiveICs,ReactionTerm = QuestionFiveSourceTerm,DiffusionConstant=DiffusionConstant)",sort="cumtime")