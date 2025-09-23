from casadi import *
from dhMat import *
from Inertia import *
class NE:
    def __init__(self,n,H,mass,g,CoM,q,dq):
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
            z_i = H[i].R() * z
            zi.append(z_i)
        for i in range(0,n):
            o_i = H[i].O();
            oi.append(o_i)
        for i in range(n):
            p_i = H[i].v * CoM[i]
            pi.append(p_i[0:3])
        

        for i in range(0,len(I)):
            I[i].RIR_T(H[i].R())
            self.M = self.M + mass[i]*(Jv[i].T@Jv[i]) +Jw[i].T@I[i].v@Jw[i]
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    cijk = 0.5*(jacobian(self.M[i,j],q[k]) + jacobian(self.M[i,k],q[j]) - jacobian(self.M[j,k],q[i]))*dq[k]
                    self.C[i,j] += cijk
                    