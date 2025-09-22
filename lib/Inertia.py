from casadi import *
class Inertia:
    def __init__(self):
        self.v = None
        self.m = None
    def box_com(self,m,x,y,z):
        I = Inertia()
        ratio = 1/12
        I.m = m
        mass = ratio * m
        ixx = mass *(y*y + z*z)
        iyy = mass *(x*x + z*z)
        izz = mass *(x*x + y*y)
        I.v = SX([[ixx, 0,   0],
                 [0,   iyy, 0],
                 [0,   0,   izz]])
        return I
    def cylinder_com(self,m,r,h):
        I = Inertia()
        I.m = m
        ratio1 = 1/12
        ratio2 = 1/2
        mass1 = ratio1*m
        mass2 = ratio2*m
        ixx = mass1*(3*r*r+h*h)
        iyy = ixx
        izz = mass2*r*r
        I.v = SX([[ixx, 0,   0],
                 [0,   iyy, 0],
                 [0,   0,   izz]])
        return I
    def parallerAxis(self,dx,dy,dz):
        if(self.v == None):
            return None
        d = SX([dx,dy,dz])
        d2 = (d.T@d)[0,0]
        self.v = self.v + self.m *(d2*SX.eye(3) - d@d.T)
    def RIR_T(self,R):
        if(self.v == None):
            return None
        self.v = R@self.v@R.T
    
