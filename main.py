from casadi import *
from lib import dhMat
from lib import Inertia
from lib import strChange

h = dhMat.dhMat()
inta = Inertia.Inertia()

q1 = SX.sym('q1')
q2 = SX.sym('q2')
q3 = SX.sym('q3')


h01 = h.RotZ(q1);
h12 = h.Trans(0, 0, 0.010)*h.RotY(pi/2) * h.RotZ(q2);
h23 = h.Trans(0.050,0,0) * h.RotZ(q3);
h3e = h.Trans(0.050, 0, 0);
h02 = h01 * h12;
h03 = h01 * h12 * h23;
h0e = h01 * h12 * h23 * h3e;
H = [h01,h02,h03]
z = SX([0,0,1.0])
e = h0e

mass1 = 3
mass2 = 5
mass3 = 5

i1 = inta.cylinder_com(mass1,0.01,0.01)
i2 = inta.box_com(mass2,0.05,0.005,0.005)
i2.parallerAxis(-0.025,0,0)
i3 = inta.box_com(mass2,0.05,0.005,0.005)
i2.parallerAxis(-0.025,0,0)
I = [i1,i2,i3]

mass = [mass1,mass2,mass3]
o0 = SX([0,0,0])
o1 = h01.v[0:3, 3]
o2 = h02.v[0:3, 3]

p1 = h01.v@SX([0,0,5,1])
p2 = h02.v@SX([25,0,0,1])
p3 = h03.v@SX([25,0,0,1])  

z0 = h01.R()@z
z1 = h02.R()@z
z2 = h03.R()@z


Jv1 = horzcat(cross(z0,p1[0:3] - o0),SX.zeros(3,1),SX.zeros(3,1))
Jv2 = horzcat(cross(z0,p2[0:3] - o0),cross(z1,p2[0:3] - o1),SX.zeros(3,1))
Jv3 = horzcat(cross(z0,p3[0:3] - o0),cross(z1,p3[0:3] - o1),cross(z2,p3[0:3] - o2))
Jv = [Jv1,Jv2,Jv3]

Jw1 = horzcat(z0, SX.zeros(3,1), SX.zeros(3,1))
Jw2 = horzcat(z0, z1, SX.zeros(3,1))
Jw3 = horzcat(z0, z1, z2)
Jw = [Jw1,Jw2,Jw3]
M = SX.zeros(3,3)
for i in range(0,len(I)):
    I[i].RIR_T(H[i].R())
    M = M + mass[i]*(Jv[i].T@Jv[i]) +Jw[i].T@I[i].v@Jw[i]
for i in range(3):
    for j in range(3):
        strM = str(M[j,i])
        strM = strChange.strChange(strM,3)
        print(strM)

