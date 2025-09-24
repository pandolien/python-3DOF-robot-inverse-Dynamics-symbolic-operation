from casadi import *
from .dhMat import *
from .Inertia import *
from .strChange import *
class NE:
    def __init__(self,n,H,mass,I,g,CoM,q,dq):
        self.M = SX.zeros(n,n)
        self.C = SX.zeros(n,n)
        self.G = SX.zeros(n,1)
        self.n = n
        z = SX([0,0,1])
        oi = []
        zi = []
        pi = []
        Jw = []
        Jv = []
        for i in range(0,n):
            z_i = H[i].R() @ z
            zi.append(z_i)
        for i in range(0,n):
            o_i = H[i].O();
            oi.append(o_i)
        for i in range(n):
            p_i = H[i].v @ CoM[i]
            pi.append(p_i[0:3])
        for i in range(n):
            Jv_i = []
            Jw_i = []
            for j in range(n):
                if(j <= i):
                    Jv_i.append(cross(zi[j],pi[i] - oi[j]))
                    Jw_i.append(zi[j])
                else:
                    Jv_i.append(SX.zeros(3,1))
                    Jw_i.append(SX.zeros(3,1))
            Jv.append(horzcat(*Jv_i))
            Jw.append(horzcat(*Jw_i))
        for i in range(0,len(I)):
            I[i].RIR_T(H[i].R())
            self.M = self.M + mass[i]*(Jv[i].T@Jv[i]) +Jw[i].T@I[i].v@Jw[i]
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    cijk = 0.5*(jacobian(self.M[i,j],q[k]) + jacobian(self.M[i,k],q[j]) - jacobian(self.M[j,k],q[i]))*dq[k]
                    self.C[i,j] += cijk
        for i in range(n):
            self.G += Jv[i].T * SX(mass[i]) * g
    def GetMassMatrix(self,i,j):
        STR = str(self.M[i,j])
        STR_ = strChange(STR)
        return STR_
    def GetCorioliMatrix(self,i,j):
        STR = str(self.C[i,j])
        STR_ = strChange(STR)
        return STR_
    def GetGravityVector(self,i):
        STR = str(self.G[i])
        STR_ = strChange(STR)
        return STR_
    
