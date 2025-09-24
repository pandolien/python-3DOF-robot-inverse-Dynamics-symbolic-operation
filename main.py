from casadi import *
from lib import dhMat
from lib import Inertia
from lib import NewtonEuler
h = dhMat.dhMat()
inta = Inertia.Inertia()

q1 = SX.sym('q1')
q2 = SX.sym('q2')
q3 = SX.sym('q3')
q = [q1,q2,q3] 

dq1 = SX.sym('dq1')
dq2 = SX.sym('dq2')
dq3 = SX.sym('dq3')

dq = [dq1,dq2,dq3]

gx = SX.sym('gx')
gy = SX.sym('gy')
gz = SX.sym('gz')
g = vertcat(gx, gy, gz)

mass1 = 3.0
mass2 = 5.0
mass3 = 5.0
mass = [mass1,mass2,mass3]

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

CoM1 = SX([0,0,0.05,1])
CoM2 = SX([0.025,0,0,1])
CoM3 = SX([0.025,0,0,1])
CoM = [CoM1,CoM2,CoM3]


i1 = inta.cylinder_com(mass1,0.01,0.01)
i2 = inta.box_com(mass2,0.05,0.005,0.005)
i2.parallerAxis(-0.025,0,0)
i3 = inta.box_com(mass2,0.05,0.005,0.005)
i2.parallerAxis(-0.025,0,0)
I = [i1,i2,i3]
ne = NewtonEuler.NE(3,H,mass,I,g,CoM,q,dq)


