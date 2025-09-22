from casadi import *
class dhMat:
    def __init__(self):
        self.v = None
    def Trans(self,x,y,z):
        h = dhMat();
        h.v = vertcat(horzcat(1,   0,  0,  x),
                        horzcat(0,   1,  0,  y),
                        horzcat(0,   0,  1,  z),
                        horzcat(0,   0,  0,  1))
        return h;
    def RotZ(self,angle):
        h = dhMat();
        h.v = vertcat(horzcat(cos(angle),   -sin(angle),  0,  0),
                        horzcat(sin(angle),   cos(angle),  0,  0),
                        horzcat(0,   0,  1,  0),
                        horzcat(0,   0,  0,  1))
        return h;

    def RotY(self,angle):
        h = dhMat();
        h.v = vertcat(horzcat(cos(angle),   0,  sin(angle),  0),
                        horzcat(0,   1,  0,  0),
                        horzcat(-sin(angle),   0,  cos(angle),  0),
                        horzcat(0,   0,  0,  1))
        return h;
    def RotX(self,angle):
        h = dhMat();
        h.v = vertcat(horzcat(1,   0,  0,  0),
                        horzcat(0,   cos(angle),  -sin(angle),  0),
                        horzcat(0,   sin(angle),  cos(angle),  0),
                        horzcat(0,   0,  0,  1))
        return h;
    def R(self):
        return self.v[0:3,0:3]
    def __mul__(self,other):
        if(self.v == None):
            return None
        h = dhMat()
        h.v = self.v@other.v
        return h
