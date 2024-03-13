#%%
import Week16General
#import Week17Functions
from Week17Functions import *

import ODESolver
Param0 = 1
FixBeta = lambda u,t :BetaFormHopf(u, t,Param0)
Ans,Time = Week16General.Shooting(FixBeta,[1.4,1],30)
print("WTF where is my answer")
print(Ans)
print(Time)

# %%
