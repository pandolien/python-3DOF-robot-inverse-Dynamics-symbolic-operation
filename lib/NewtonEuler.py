from casadi import *
from dhMat import *
from Inertia import *
class NE:
    def __init__(self,H,Jv,Jw,I,mass):
        self.M = SX.zeros(3,3)
        self.C = 
        for i in range(0,len(I)):
            I[i].RIR_T(H[i].R())
            self.M = self.M + mass[i]*(Jv[i].T@Jv[i]) +Jw[i].T@I[i].v@Jw[i]