from scipy.optimize import root
import numpy as np
def printer(x0,q):
    print(x0)
    print(q)
    return
def TestingFunc(x,y,z):
    printer(x,y)
    def shouter(x,q=y):
        print(q)
        print(x,x)
        
        return 0
    shouter(x)
    return np.array([x[0]-y**2-z**2,x[1]],dtype=np.float32)

#y=1
#z=1
A = TestingFunc([2,2],1,1)#root(TestingFunc,[1,2], args=(2,2)).x
shouter(2)
print(A)
