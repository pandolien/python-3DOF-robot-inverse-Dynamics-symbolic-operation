from casadi import *
from dhMat import *
from Inertia import *
class NE:
    def __init__(self,n,H,Jv,Jw,I,mass,g,q,dq,CoM):
        self.M = SX.zeros(n,n)
        self.C = SX.zeros(n,n)
        self.G = SX.zeros(n,1)
        self.n = n
        for i in range(0,len(I)):
            I[i].RIR_T(H[i].R())
            self.M = self.M + mass[i]*(Jv[i].T@Jv[i]) +Jw[i].T@I[i].v@Jw[i]
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    cijk = 0.5*(jacobian(self.M[i,j],q[k]) + jacobian(self.M[i,k],q[j]) - jacobian(self.M[j,k],q[i]))*dq[k]
                    self.C[i,j] += cijk
                    