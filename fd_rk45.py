from numpy import array
from math import *
from numpy import linalg
from numpy import arange
import pygame
from pygame.locals import *
#
# ----------------------------------------
#
def declare_matrix(n, m):
    #
    return array([[0.0 for j in range(0, m)] for i in range(0, n)]);


#
# ----------------------------------------
#
def declare_vector(n):
    #
    return array([0.0 for i in range(0, n)]);


#
# ----------------------------------------
#
def Jacob_AB(RHS, y, t, u_control, n, m):
    #
    A = declare_matrix(n, n);
    B = declare_matrix(n, m);
    dy = 1.0e-6;
    f0 = RHS(y, t, u_control);
    for i in range(0, n):
        yp = array(y);
        yp[i] += dy;
        f = RHS(yp, t, u_control);
        for j in range(0, n):
            A[j, i] = (f[j] - f0[j]) / dy;

    for i in range(0, m):
        up = array(u_control);
        up[i] += dy;
        f = RHS(y, t, up);
        for j in range(0, n):
            B[j, i] = (f[j] - f0[j]) / dy;

    return A, B;


#
# ----------------------------------------
#
def RHS(x, t, u_control):
    #
    dx_0 = x[1];
    # dx_1 = -sin(x[0]) + u_control[0];
    dx_1 = -x[0] + u_control[0];
    x_prim = array([dx_0, dx_1]);

    return x_prim;


#
#
# ----------------------------------------
#
# ----------------------------------------
def rhs(x, t, u):
    n=len(x)
    dx_dt=declare_vector(n)

    az_turbulence = 0.0;
    ax_wind = 0.0;
    az_wind = 0.0;

    deg2rad = pi / 180.0;
    g = 9.81;

    S = 1.0;

    mass = 25.0;
    Iy = 100;

    vx = x[0] + ax_wind;
    vz = x[1] + az_wind;

    alpha = atan2(vz, vx);
    V = sqrt(vz * vz + vx * vx);

    CD_0 = 0.30;
    CD = CD_0;

    ro_0 = 1.225;
    ro = ro_0 * (1.0 - abs(x[4]) / 44300.0) ** 4.256;

    Q_dyn = 0.5 * ro * V * V;

    L = 0.0;
    D = Q_dyn * S * CD;
    G = mass * g;

    Th = 1.0;

    Thrust_1 = 0.5 * G + u[0];
    Thrust_2 = 0.5 * G + u[1];


    Thrust_1 = x[6]
    Thrust_2 = x[7]

    cm_q = -0.01;

    Tau = 0.05;

    beta = 0.0 * deg2rad;
    cb = cos(beta);
    sb = sin(beta);

    dx_dt[0] = (-D * cos(alpha) + L * sin(alpha) - G * sin(x[5]) - Thrust_1 * sb + Thrust_2 * sb) / mass - x[2] * vz;
    dx_dt[1] = (-D * sin(alpha) - L * cos(alpha) + G * cos(x[5]) - Thrust_1 * cb - Thrust_2 * cb) / mass + x[2] * vx + az_turbulence;
    dx_dt[2] = (0.5 * (Thrust_2 * cb - Thrust_1 * cb) + cm_q * x[2]) / Iy;
    dx_dt[3] = cos(x[5]) * vx + sin(x[5]) * vz;
    dx_dt[4] = -sin(x[5]) * vx + cos(x[5]) * vz;
    dx_dt[5] = x[2];
    dx_dt[6] = (1.0/Tau)*( -x[6] + Th*u[0] );
    dx_dt[7] = (1.0/Tau)*( -x[7] + Th*u[1] );

    return dx_dt

#
#
# ----------------------------------------


def fd_rk45( RHS , x , t , dt , u_control ):
#
  y0 = RHS( x , t, u_control );

  t1 = t + dt*0.25;
  vec = x + dt*(0.25)*y0;
  y1 = RHS( vec , t1, u_control );

  t2 = t + dt*(3.0/8.0);
  vec = x + dt*( (3.0/32.0)*y0 + (9.0/32.0)*y1 );
  y2 = RHS( vec , t2, u_control );

  t3 = t + dt*(12.0/13.0);
  vec = x + dt*( (1932.0/2197.0)*y0 + (-7200.0/2197.0)*y1 + (7296.0/2197.0)*y2 );
  y3 = RHS( vec , t3, u_control );

  t4 = t + dt;
  vec = x + dt*( (439.0/216.0)*y0 + (-8.0)*y1 + (3680.0/513.0)*y2 + (-845.0/4104.0)*y3 );
  y4 = RHS( vec , t4, u_control );

  t5 = t + dt*(1.0/2.0);
  vec = x + dt*( -(8.0/27.0)*y0 + (2.0)*y1 + (-3544.0/2565.0)*y2 + (1859.0/4104.0)*y3 + (-11.0/40.0)*y4 );
  y5 = RHS( vec , t5, u_control );

  y = x + dt * ( (16.0/135.0)*y0 + (6656.0/12825.0)*y2 + (28561.0/56430.0)*y3 + (-9.0/50.0)*y4 + (2.0/55.0)*y5 );

  return y
#


def fd_rk45a(x , t , dt , u ):
#

  y=declare_vector(len(x))
  y0 = rhs(x,t,u);

  t1 = t + dt*0.25;
  vec = x + dt*(0.25)*y0;
  y1 = rhs( vec , t1, u);

  t2 = t + dt*(3.0/8.0);
  vec = x + dt*( (3.0/32.0)*y0 + (9.0/32.0)*y1 );
  y2 = rhs( vec , t2, u );

  t3 = t + dt*(12.0/13.0);
  vec = x + dt*( (1932.0/2197.0)*y0 + (-7200.0/2197.0)*y1 + (7296.0/2197.0)*y2 );
  y3 = rhs( vec , t3, u );

  t4 = t + dt;
  vec = x + dt*( (439.0/216.0)*y0 + (-8.0)*y1 + (3680.0/513.0)*y2 + (-845.0/4104.0)*y3 );
  y4 = rhs( vec , t4, u );

  t5 = t + dt*(1.0/2.0);
  vec = x + dt*( -(8.0/27.0)*y0 + (2.0)*y1 + (-3544.0/2565.0)*y2 + (1859.0/4104.0)*y3 + (-11.0/40.0)*y4 );
  y5 = rhs( vec , t5, u );

  y = x + dt * ( (16.0/135.0)*y0 + (6656.0/12825.0)*y2 + (28561.0/56430.0)*y3 + (-9.0/50.0)*y4 + (2.0/55.0)*y5 );
  #print("y=",y)
  return y
#



#

def aa_trajectory(X,Vel,dt):
    Z = 1.0;
    ###!! X to  droga do przodu ten trapez to pokazuje jak wygląda teren
    dx=0.1

    if (X <= 1.0):
        Z = 1.0;

    if (X > 1.0 and X < 1.5):
        #Z = 1.0 + (X - 1.0) * 10.0;
        Z=0.5
    if (X >= 1.5 and X <= 2.0):
        #Z = 6.0;
        Z = 6.0 - (X - 2.0) * 10.0;
        Z=1
    if (X >= 2.0 and X <= 2.5):
        Z = 6.0 - (X - 2.0) * 10.0;
        #Z=1
    if (X >= 2.5):
        Z = 1.0;


    Z1 = Z;

    dx = Vel * dt;
    X = X + dx;


    if (X <= 1.0):
        Z = 1.0;

    if (X > 1.0 and X < 1.5):
        Z = 1.0 + (X - 1.0) * 10.0;
        Z=0.5
    if (X >= 1.5 and X <= 2.0):
        Z = 6.0;
        Z=1
    if (X >= 2.0 and X <= 2.5):
        Z = 6.0 - (X - 2.0) * 10.0;

    if (X >= 2.5):
        Z = 1.0;

    alpha = atan2(Z - Z1, dx);

    return Z,alpha

#
# ----------------------------------------
#
def aa_matrices_AB(x,t,u,n, m):
    de = 0.00001
    f0= rhs( x , t, u )
    #print("f0=",f0)

    A=declare_matrix(n,n)
    B=declare_matrix(n,m)
    #print("ssssssdadsad")
    for j in range(0,n):
        dx=declare_vector(n)
        dx[j]=de
        #print("ffffffffffffffffffffffff")
        A[:,j]=(rhs(x+dx,t,u)-f0)/de
        #print("ffffffffffffffffffffffff")
        #print(j)
    for j in range(0,m):
        du=declare_vector(m)


        du[j]=de
        B[:,j]=(rhs(x,t,u+du)-f0)/de
    #print("ssssssdadsad")
    return A, B

#
# ----------------------------------------
#
def aa_euler( RHS , x , t , dt, u ):

    c_bet = 0.90;
    del_state = 1.0e-6;
    coef = - c_bet / del_state;
    n = len(x);
    y=declare_vector(n)
    Jacob=declare_matrix(n,n)

    y0 = rhs(x, t, u);
    for i in range(0,n):
        tmp = x[i];
        x[i] += del_state;
        vec = rhs(x, t, u);
        x[i] = tmp;
        Jacob[:, i] = coef * (vec - y0);

    on=declare_matrix(n,n)

    for i in range(0,n):

        on[i,i]=1

    Jacob += on/ dt;
    dx = linalg.inv(Jacob).dot(y0);
    y = x + dx;
    #print("on", y)
    return y


#
# ----------------------------------------
#
def declare_trajectory():

    YELLOW = (255, 255, 0)
    BLUE = (0, 0, 0)

    green = (0, 0, 0)
    blue = (255, 255, 255)


    pygame.init()

    # Set dimensions of game GUI
    w, h = 640, 750    #!!!!!!!!!!!!!!!!!!!!! jakto zmieniasz to też w tej  drugiej funckji transfo.... coś tam
    screen = pygame.display.set_mode((w, h))

    y = []
    x = []
    pole = []
    pop_s=[]
    G=[0,0]
    H = []
    F = []
    yd = []
    xd=[]
    # Set running and moving values
    running = True

    k = (100.5, 100.5)
    # Setting what happens when game
    # is in running state
    screen.fill(YELLOW)
    koniec = 0

    d = 10
    pop = (0, 500)  # poprzedni srodek/miejsce startowe
    # pygame.draw.line(screen, BLUE, (40, 0), (40, 100), 1)
    pygame.draw.rect(screen, (200, 0, 0), (pop, (d, d)), 5)
    for i in range(0, w, d):
        for j in range(0, h, d):
            r = pygame.draw.rect(screen,YELLOW, (i, j, d, d), 1)  #kratownica

            pole.append(r)
        ###################################################################kratownica

        for i in range(0,h,100):
            pygame.draw.line(screen, BLUE, (0, h-i), (w,h-i), 1)
            pygame.draw.line(screen, BLUE, (i, 0), (i, h), 1)

        #######################################################

        font = pygame.font.Font('freesansbold.ttf', 20)

        # create a text surface object,
        # on which text is drawn on it.

        text2 = font.render("0", True, green, blue)
        textRect2 = text2.get_rect()
        textRect2.center = (8, 740)
        screen.blit(text2, textRect2)

        text2 = font.render("100", True, green, blue)
        textRect2 = text2.get_rect()
        textRect2.center = (100, 740)
        screen.blit(text2, textRect2)

        text2 = font.render("200", True, green, blue)
        textRect2 = text2.get_rect()
        textRect2.center = (200, 740)
        screen.blit(text2, textRect2)

        text2 = font.render("300", True, green, blue)
        textRect2 = text2.get_rect()
        textRect2.center = (300, 740)
        screen.blit(text2, textRect2)

        text2 = font.render("400", True, green, blue)
        textRect2 = text2.get_rect()
        textRect2.center = (400, 740)
        screen.blit(text2, textRect2)

        text2 = font.render("500", True, green, blue)
        textRect2 = text2.get_rect()
        textRect2.center = (500, 740)
        screen.blit(text2, textRect2)

        text2 = font.render("600", True, green, blue)
        textRect2 = text2.get_rect()
        textRect2.center = (600, 740)
        screen.blit(text2, textRect2)
        ##################################pion
        text2 = font.render("100", True, green, blue)
        textRect2 = text2.get_rect()
        textRect2.center = (17, h-100)
        screen.blit(text2, textRect2)

        text2 = font.render("200", True, green, blue)
        textRect2 = text2.get_rect()
        textRect2.center = (17, h-200)
        screen.blit(text2, textRect2)

        text2 = font.render("300", True, green, blue)
        textRect2 = text2.get_rect()
        textRect2.center = (17, h-300)
        screen.blit(text2, textRect2)

        text2 = font.render("400", True, green, blue)
        textRect2 = text2.get_rect()
        textRect2.center = (17, h-400)
        screen.blit(text2, textRect2)

        text2 = font.render("500", True, green, blue)
        textRect2 = text2.get_rect()
        textRect2.center = (17, h-500)
        screen.blit(text2, textRect2)

        text2 = font.render("600", True, green, blue)
        textRect2 = text2.get_rect()
        textRect2.center = (17, h-600)
        screen.blit(text2, textRect2)

        text2 = font.render("700", True, green, blue)
        textRect2 = text2.get_rect()
        textRect2.center = (17, h - 700)
        screen.blit(text2, textRect2)
        pygame.display.set_caption('Narysuj kontury krajobraz wybierając odpowiednie punkty ')
        #######################################################




    while running:

        for event in pygame.event.get():

            # Close if the user quits the
            # game

            if event.type == QUIT:
                running = False

            # Making the image move
            elif event.type == MOUSEBUTTONDOWN:
                # if rect.collidepoint(event.pos):
                for i in range(1, len(pole)):
                    if pole[i - 1].collidepoint(event.pos):
                        if pop != pole[i - 1].center:

                            pygame.draw.rect(screen, BLUE, pole[i - 1], 4)
                            col = pygame.draw.line(screen, BLUE, pop, pole[i - 1].center, 6)
                            if abs(pop[0]) >= abs(pole[i - 1].center[0]):
                                running = False

                            for j in range(abs(pop[0]), abs(pole[i - 1].center[0])):

                                y = (h-pole[i - 1].center[1] - h-pop[1]) * (j - pop[0]) / (pole[i - 1].center[0] - pop[0]) + \
                                    h-pop[1]
                                yd.append(y)
                                xd.append(j)

                            pop_s.append(pop)

                            pop = pole[i - 1].center


                k = event.pos


        pygame.display.update()

    pygame.quit()

    return pop_s,xd;


#
# ----------------------------------------
#
def aa_trajectory_transformation(X,xd,pop,Vel,dt):
    #Z = 1.0;
    ###!! X to  droga do przodu ten trapez to pokazuje jak wygląda teren
    w, h = 640, 750  # !!!!!!!!!!!!!!!!!!!!! jakto zmieniasz to też w tej  drugiej funckji transfo.... coś tam
    dx=0.1
    K=1
    ii=0
    KD=1
    koniec=0

    for i in range(0,len(pop)-1):
        if pop[i][0] <= X*K and X*K<pop[i+1][0]:
            ii=i

        elif X*K >=pop[len(pop)-1][0]:

            koniec=1

    Z = (pop[ii + 1][1] - pop[ii][1]) * (X*KD- pop[ii][0]) / (pop[ii + 1][0] - pop[ii][0]) + pop[ii][1]

    Z1 = h-Z;

    dx = Vel * dt;
    X = X + dx;

    for i in range(0,len(pop)-1):
        if pop[i][0] <= X*K and X*K<pop[i+1][0]:
            ii=i
        elif X * K >= pop[len(pop)-1][0]:

            koniec = 1
    Z = (pop[ii + 1][1] - pop[ii][1]) * (X*KD - pop[ii][0]) / (pop[ii + 1][0] - pop[ii][0]) + pop[ii][1]



    alpha = atan2(h-Z - Z1, dx);

    return h-Z,alpha,koniec


