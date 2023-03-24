from control import lqr
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
import pygame
from pygame.locals import *
from fd_rk45 import *

global K;

global az_turbulence;
global ax_wind;
global az_wind;


tp = [];
xp0 = [];
xp1 = [];
up = [];
gp = [];
zp = [];

plt.figure(figsize=(15, 7))
plt.ion()
ax = plt.gca()


line0, = ax.plot([], [], 'r')
line1, = ax.plot([], [], 'b')
lineu, = ax.plot([], [], '-y')


handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels)

n = 8;
m = 2;

z0 = 250.0;
h_flight = 10.0;  # Over the terrain
c_turb = 1000.0;
X_turb_1 = 1500.0;
X_turb_2 = 2000.0;

t_end = 352.02;
#t_end=30

x = declare_vector(n);
#x = [10, 10];
x = [ 0.0, 0.0, 0.0, 0.0, -z0, 0.0, 0.0, 0.0] #x=[vx,vz,wy,X,Z,kot,silnik 1,silnik 2)



Vel = 3; # / 3.6;

dt = 0.01;



xref = declare_vector(n);
xref = [1, 0];
u = declare_vector(m);
u = [ 0.0 , 0.0 ];

e=declare_vector(n)
#print(e)
t = 0;
#t=0

#Q = array([[100000,0],[0,1]])
R = array([[1,0],[0,1]])
Q=declare_matrix(n,n)
"""
Q[0,0]=1000
Q[1,1]=100000
Q[2,2]=1000
Q[3,3]=1
Q[4,4]=1000
Q[5,5]=100000
Q[6,6]=10000
Q[7,7]=10000
Q=Q*0.01
"""
Q[0,0]=10000000
Q[1,1]=10000000
Q[2,2]=1000
Q[3,3]=100000
Q[4,4]=1000000
Q[5,5]=10000000
Q[6,6]=10000
Q[7,7]=10000
Q=Q*1

R=R*1
r=[0.1]
#print("Q=",Q)

#A, B = Jacob_AB(RHS, x, t, u_control, n, m);

#A=[[0,1],[-9.81,0]]
#B=[[0],[1]]

#A=[[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[20,0],[0,20]]
#B=[[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[20,0],[0,20]]

#y2 = rhs(x, t, u)

################################### window

pop_s,xdw=declare_trajectory()
print("yd",pop_s)
print("yd",pop_s[1])
print("yd",pop_s[1][0])
#print("xd",xdw)
###################################

while t <= t_end + dt:
    #A=4/0
    A, B = aa_matrices_AB(x, t, u, n, m)
    # B[7,1]=20
    """
    for i in range(0,n-1):
        A[i,1]=20

    """
    #print("u=",u)
    #print("A=",A)
    #print("sasasadsasasasasasasaassa")
    #K = lqr(A, B, Q, R)[0];

    X=x[3]
    Z0=5
    #z_terr,alfa=aa_trajectory(X,Vel,dt)
    z_terr,alfa,koniec=aa_trajectory_transformation(X, xdw, pop_s, Vel, dt)
    if koniec==1:
        break;
    Vx = Vel * cos(alfa);
    Vz = Vel * sin(alfa);

    z_ref = z_terr + h_flight;
    x_ref=X
    #R=[]
   # A, B = Jacob_AB(RHS, x, t, u_control, n, m);

    e = declare_vector(n)
    """
    e[1] = x[1] - (cos(x[6]) * Vx + sin(x[6]) * Vz)
    e[2] = x[2] - (sin(x[6]) * Vx - cos(x[6]) * Vz)
    e[3] = x[3] - 0;
    e[4] = x[4] - (x_ref);
    e[4] = 0.0;
    e[5] = x[5] - (-z_ref);
    e[6] = x[6] - 0.0;
    """
    e[0] = x[0] - (cos(x[5]) * Vx + sin(x[5]) * Vz)
    e[1] = x[1] - (sin(x[5]) * Vx - cos(x[5]) * Vz)
    e[2] = x[2] - 0;
    e[3] = x[3] - (x_ref);
    e[3] = 0.0;
    e[4] = x[4] - (-z_ref);
    e[5] = x[5] - 0.0;

    #print("e=",e)
    #print("K=",K)

    K = lqr(A, B, Q, R)[0];
    #e = [x[0]-xref[0], x[1]-xref[1]];  #e=x-xref

    u = -K.dot(e);

    ######po co te uszukiwanie
    umax=1000000
    u[0]=max(-umax,min(umax,u[0]))
    u[1]=max(-umax,min(umax,u[1]))
    #########################
    #print("u=",u)

    tp.append(X);
    xp0.append(x[0]);     #yp
    xp1.append(x[1]);       #yp
    up.append(z_terr);
    ######
    gp.append(z_ref);
    zp.append(-x[4]);
    """
    #X[2] W , X[6]T1,X[7]T2
    tp.append(t);
    xp0.append(x[0]);  # yp
    xp1.append(x[1]);  # yp
    up.append(x[6]);
    ######
    gp.append(x[7]);
    zp.append(1000*x[2]);
    """

    #txt = 't = {:8.4f}   x = {:9f}   v = {:9f}   u = {:9f}'.format(t, x[0], x[1], u[0]);
    #line0.set_data(tp, zp);
    #line1.set_data(tp, gp);                        #line1.set_data(tp, xp1);
    #lineu.set_data(tp, up);
    line0.set_data(tp, zp);
    line1.set_data(tp, gp);  # line1.set_data(tp, xp1);
    lineu.set_data(tp, up);
    #print("=",gp)
    #print("X=", X)
    """
    if t > 50.0:
        tp.pop(0);
        xp0.pop(0);
        xp1.pop(0);
        up.pop(0);
    """
    ##x = fd_rk45(RHS, x, t, dt, u);
    #x = fd_rk45a(x, t, dt, u);
    x=aa_euler(RHS, x, t, dt, u)
    print("t=",t)
    #print("sssdnowyyyyyyyyyyyyyyyyyyyyyyyyyyyyyksajdhkssss")
    t += dt;

ax.relim()
ax.autoscale_view(False, True, False)
ax.grid(True);
legend([line0, line1, lineu], ["trajektoria drona", "trajektoria referencyjną", "ukształtowanie ziemi"], loc=2);
plt.subplots_adjust(left=0.07, right=0.95, bottom=0.1, top=0.87);
plt.title("Calkowanie rownania wahadla\n\n");
plt.draw();

plt.pause(199999993.0)


