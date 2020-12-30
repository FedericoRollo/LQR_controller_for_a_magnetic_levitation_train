%% WHOLE SIMULATION

clear all
clc


%% ----------------------- HEADWAY motion -------------------------------

z=0.05; % nominal z

%motor constants
i_e=100; % nominal current which provides nominal z (it MUST be replaced with the REAL NOMINAL VALUE)
Rx=1.04; %Rxesistance of the stator
Ld=0.06; %inductance
Lq=0.06; %inductance
M=0.1; %mutual inductance stator rotor
tau=1; %pole pitch
p=3;% poles

%train constants
m=60000;%mass 
rho=1.225; %aerodynamic coefficient
Sa=13; %aerodynamic coefficient
Cd=1.2; %aerodynamic coefficient


id1_ref=0;
id2_ref=0;

xdot_ref=430/3.6; %m/s


% stators currents

syms iq1 iq2 id1 id2 xdot Vq1 Vq2 Vd1 Vd2

states = [iq1; id1; iq2; id2; xdot];
inputs = [Vq1; Vd1; Vq2; Vd2];

iq1_dot = -(Rx/Lq)*iq1 + (Ld/Lq)*((M/Ld)*i_e + id1)*xdot*pi/tau + Vq1/Lq ;
iq2_dot = -(Rx/Lq)*iq2 + (Ld/Lq)*((M/Ld)*i_e + id2)*xdot*pi/tau + Vq2/Lq ;

id1_dot = -(Rx/Ld)*id1 - (Lq/Ld)*(iq1)*xdot*pi/tau + Vd1/Ld ;
id2_dot = -(Rx/Ld)*id2 - (Lq/Ld)*(iq2)*xdot*pi/tau + Vd2/Ld ;

Fx1 = (3/2)*(pi/tau)*(Ld*((M/Ld)*i_e) + (Ld-Lq)*id1)*iq1;
Fx2 = (3/2)*(pi/tau)*(Ld*((M/Ld)*i_e) + (Ld-Lq)*id2)*iq2;

xdot_dot = (1/m)*(p*Fx1 + p*Fx2 - (rho*Sa*Cd/2)*xdot^2); % + vy*wpsi - vz*wtheta

ode = [iq1_dot; id1_dot; iq2_dot; id2_dot; xdot_dot];



%references for currents, related to speed
Fx1_ref = ((rho*Sa*Cd/2)*xdot_ref^2)/(2*p);
Fx2_ref = Fx1_ref;

iq1_ref = (Fx1_ref)/((3/2)*(pi/tau)*(Ld*(M/Ld)*i_e));
iq2_ref = (Fx2_ref)/((3/2)*(pi/tau)*(Ld*(M/Ld)*i_e));

Vq1_ref = Rx*iq1_ref-Ld*((M/Ld)*i_e+id1_ref)*xdot_ref*(pi/tau);
Vq2_ref = Rx*iq2_ref-Ld*((M/Ld)*i_e+id2_ref)*xdot_ref*(pi/tau);

Vd1_ref = Rx*id1_ref+Lq*(iq1_ref)*xdot_ref*pi/tau ;
Vd2_ref = Rx*id2_ref+Lq*(iq2_ref)*xdot_ref*pi/tau ;

Ax=jacobian(ode, states);
Ax=subs(Ax, {iq1 id1 iq2 id2 xdot}, {iq1_ref id1_ref iq2_ref id2_ref xdot_ref});
Ax=eval(Ax);

Bx= jacobian(ode, inputs);
Bx=subs(Bx, {Vq1 Vd1 Vq2 Vd2}, {Vq1_ref Vd1_ref Vq2_ref Vd2_ref});
Bx=eval(Bx);

Cx = [1 0 0 0 0;
      0 1 0 0 0;
      0 0 0 1 0];

Dx = [0 0 0 0;
      0 0 0 0; 
      0 0 0 0];

Qx=[1000 0 0 0 0; 0 100000 0 0 0; 0 0 1000 0 0; 0 0 0 100000 0; 0 0 0 0 10];
Rrx=[10 0 0 0; 0 10 0 0; 0 0 10 0; 0 0 0 10];
Rrkx=10*eye(3);

Xref = [iq1_ref; id1_ref; iq2_ref; id2_ref; xdot_ref];
Vref_x = [Vq1_ref; Vd1_ref; Vq2_ref; Vd2_ref];

X_init = [0 ;0;0;0;0];

[Kx,Sx,ex]=lqr(Ax,Bx,Qx,Rrx);

%% ------------------- SWAY MOTION ------------------------------

% Constants

radius = 0;
vx=430/3.6;
xB=8;
yB=2;
mu0=1.25*10^-6;
Ng=500;
Ag=1;
m=60000;
g=9.8;
Ra=0;
Rg=100;
Jz=2853*10^3;
phiy=0;
wtheta=0;
P1=100;
P2=100;
P3=100;
P4=100;
P5=100;
P6=100;
ymax = 0.02;
l=8;

[shape,int11,int21,int31,omega_1]=Shape_Function();
w=eval(omega_1);       %frequenza naturale modo vibrazionale(provvisorio)
PHI=shape;     %shape function?
PHIf=eval(int11);    %moltiplicativi forze per vibrazioni
PHIc=eval(int21);
PHIr=eval(int31);


% ------------------------STATES STATEMENT--------------------------------
syms y vy psi wpsi ql vql ig1 ig2 ig3 ig4 ig5 ig6 

states = [ y vy psi wpsi ig1 ig2 ig3 ig4 ig5 ig6 ql vql ];
% states = [ y vy psi wpsi ig1 ig2 ig3 ig4 ig5 ig6 ];

% -----------------------INPUT STATEMENT----------------------------------
syms Vg1 Vg2 Vg3 Vg4 Vg5 Vg6 
input = [Vg1 Vg2 Vg3 Vg4 Vg5 Vg6];

s=4;
PHI4=eval(PHI);
s=12;
PHI12=eval(PHI);
s=20;
PHI20=eval(PHI);

deltay1= ymax - y - psi*xB - ql*PHI4;
deltay2= ymax - y - ql*PHI12;
deltay3= ymax - y + psi*xB - ql*PHI20;
deltay4= ymax + y + psi*xB + ql*PHI4;
deltay5= ymax + y + ql*PHI12;
deltay6= ymax + y - psi*xB + ql*PHI20;

deltaydot1= - vy - wpsi*xB - vql*PHI4;
deltaydot2= - vy - vql*PHI12;
deltaydot3= - vy + wpsi*xB - vql*PHI20;
deltaydot4= - deltaydot1; 
deltaydot5= - deltaydot2;
deltaydot6= - deltaydot3;

% guidance forces
G1 = (mu0*Ng^2*Ag/4)*(ig1/deltay1)^2;
G2 = (mu0*Ng^2*Ag/4)*(ig2/deltay2)^2;
G3 = (mu0*Ng^2*Ag/4)*(ig3/deltay3)^2;
G4 = (mu0*Ng^2*Ag/4)*(ig4/deltay4)^2;
G5 = (mu0*Ng^2*Ag/4)*(ig5/deltay5)^2;
G6 = (mu0*Ng^2*Ag/4)*(ig6/deltay6)^2;

% distribuzioni delle forze sul modello elastico
gf=(G4-G1 + (m*g*sin(phiy)-m*radius*cos(phiy)*wpsi^2)/3)/l;
gc=(G5-G2 +(m*g*sin(phiy)-m*radius*cos(phiy)*wpsi^2)/3)/l;
gr=(G6-G3 +(m*g*sin(phiy)-m*radius*cos(phiy)*wpsi^2)/3)/l;

%----------differential equations-----------
psidot  = wpsi;
wpsidot = ((xB*(-G4+G1+G6-G3)+yB*(-P1+P2-P3+P4-P5+P6))/Jz); %(Jy-Jz)*tethadot*psidot -> si toglie
ydot  = vy;
vydot=((-G4+G1-G5+G2-G6+G3+m*g*sin(phiy)-m*radius*(wpsi^2)*cos(phiy))/m) - vx*psidot;% + vz*phidot -> si toglie

%equazioni elastiche
qldot=vql;
vqldot=-ql*w^2+gf*PHIf+gc*PHIc+gr*PHIr;

%LKT DEGLI ELETTROAMGNETI DELLA GUIDANCE
igdot1=ig1*deltaydot1/deltay1-(2/(mu0*(Ng^2)*Ag))*deltay1*(Rg*ig1-Vg1);
igdot2=ig2*deltaydot2/deltay2-(2/(mu0*(Ng^2)*Ag))*deltay2*(Rg*ig2-Vg2);
igdot3=ig3*deltaydot3/deltay3-(2/(mu0*(Ng^2)*Ag))*deltay3*(Rg*ig3-Vg3);
igdot4=ig4*deltaydot4/deltay4-(2/(mu0*(Ng^2)*Ag))*deltay4*(Rg*ig4-Vg4);
igdot5=ig5*deltaydot5/deltay5-(2/(mu0*(Ng^2)*Ag))*deltay5*(Rg*ig5-Vg5);
igdot6=ig6*deltaydot6/deltay6-(2/(mu0*(Ng^2)*Ag))*deltay6*(Rg*ig6-Vg6);


%diffetential equation
ode = [ydot; vydot; psidot; wpsidot; igdot1; igdot2; igdot3; igdot4; igdot5; igdot6; qldot; vqldot];

% Linearization

%compute the jacobian using the symbolic variables
a = jacobian(ode,states);
b = jacobian(ode,input);

% references
y_ref  = 0;
vy_ref     = 0;
psi_ref    = 0;
wpsi_ref   = 0;
ql_ref     = 0;
vql_ref    = 0;
ig_ref = 10;

% % e.p. for linearization
y      = 0;
vy     = 0;
psi    = 0;
wpsi   = 0;
ql     = 0;
vql    = 0;
ig1    = ig_ref;
ig2    = ig_ref;
ig3    = ig_ref;
ig4    = ig_ref;
ig5    = ig_ref;
ig6    = ig_ref;

%set inputs linearization values
Vg_ref = Rg*ig_ref;
Vg1 = Vg_ref;
Vg2 = Vg_ref;
Vg3 = Vg_ref;
Vg4 = Vg_ref;
Vg5 = Vg_ref;
Vg6 = Vg_ref;

Vref_y = [ Vg1; Vg2; Vg3; Vg4; Vg5; Vg6];

%compute A B matricies
Ay=eval(a);
By=eval(b);

Cy = [1 0 0 0 0 0 0 0 0 0 0 0;
      0 0 1 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 1 0 0 0 0 0;
      0 0 0 0 0 0 0 1 0 0 0 0;
      0 0 0 0 0 0 0 0 1 0 0 0;
      0 0 0 0 0 0 0 0 0 1 0 0;
      0 0 0 0 0 0 0 0 0 0 1 0;
      0 0 0 0 0 0 0 0 0 0 0 1];
  
  Dy =[0 0 0 0 0 0; 
      0 0 0 0 0 0; 
      0 0 0 0 0 0;
      0 0 0 0 0 0; 
      0 0 0 0 0 0;
      0 0 0 0 0 0;
      0 0 0 0 0 0;
      0 0 0 0 0 0];
  
contry = ctrb(Ay,By); % controllability matrix
ry = rank(contry);

Qy = diag([1000 0.1 10 10 100000 100000 100000 100000 100000 100000 0.1 0.1]);
Rry = 1*eye(6);
Rrky=1*eye(8);

Ky = lqr(Ay,By,Qy,Rry);

%INITIAL VALUES

y_init      = 0.01;
vy_init     = 0;
psi_init    = 0;
wpsi_init   = 0;
ql_init     = 0;
vql_init    = 0;
ig1_init    = 0;
ig2_init    = 0;
ig3_init    = 0;
ig4_init    = 0;
ig5_init    = 0;
ig6_init    = 0;


Yref = [y_ref; vy_ref; psi_ref; wpsi_ref; ig_ref; ig_ref; ig_ref; ig_ref; ig_ref; ig_ref; ql_ref; vql_ref];

Y_init = [y_init;0;0;0;0;0;0;0;0;0;0;0];
%% Pump motion
% vertical motion (along z axes) is here analyzed. current 1 -> right front; 2 -> right center; 3-> right rear; 4 -> left front; 5 -> left center; 6 -> left rear   

% Constants
h = 1; %constant, height of z at rest
zi_rip = 0.05; %zi constant of zi when the train is at rest: graphically z(h)|zi(-zi_rest), where |->axes
xdot_ref = 430/3.6; % headway speed
m = 60000;
radius = 0;
xB=8;
yB=2;
zB=-0.5;
mu0=1.25*10^-6;
Nr=300;
Al=1;
rho=1.225;
Ss=13;
Cd=1.2;
Cl=0.3;
g=9.81;
Rz=100;
xA=12;
zA=0.5;
P1=100; 
P2=100;
P3=100; 
P4=100; 
P5=100; 
P6=100;
G1=10;
G2=10;
G3=10;
G4=10;
G5=10;
G6=10;
vx=400/3.6;
Jx=296*10^3;
Jy=2853*10^3;
l=8;  

[shape,int11,int21,int31,omega_1]=Shape_Function();
w=eval(omega_1);       %frequenza naturale modo vibrazionale(provvisorio)
PHI=shape;     %shape function?
PHIf=eval(int11);    %moltiplicativi forze per vibrazioni
PHIc=eval(int21);
PHIr=eval(int31);

s=4;
PHI4=eval(PHI);
s=12;
PHI12=eval(PHI);
s=20;
PHI20=eval(PHI);

% Differential equation
 
syms z zdot phi phidot theta thetadot ir1 ir2 ir3 ir4 ir5 ir6 qv vqv Vr1 Vr2 Vr3 Vr4 Vr5 Vr6 real

states = [z; zdot; phi; phidot; theta; thetadot; ir1; ir2; ir3; ir4; ir5; ir6; qv; vqv];
inputs = [Vr1; Vr2; Vr3; Vr4; Vr5; Vr6];

% vertical displacements, measured by sensors
zFL = h + zi_rip -z - phi*yB + theta*xB - qv*PHI4;
zCL = h + zi_rip -z - phi*yB - qv*PHI12;
zRL = h + zi_rip -z - phi*yB - theta*xB - qv*PHI20;
zFR = h + zi_rip -z + phi*yB + theta*xB - qv*PHI4;
zCR = h + zi_rip -z + phi*yB - qv*PHI12;
zRR = h + zi_rip -z + phi*yB - theta*xB - qv*PHI20;

% vertical speed
zFL_dot = -zdot - phidot*yB + thetadot*xB - vqv*PHI4;
zCL_dot = -zdot - phidot*yB - vqv*PHI12;
zRL_dot = -zdot - phidot*yB - thetadot*xB - vqv*PHI20;
zFR_dot = -zdot + phidot*yB + thetadot*xB - vqv*PHI4;
zCR_dot = -zdot + phidot*yB - vqv*PHI12;
zRR_dot = -zdot + phidot*yB - thetadot*xB - vqv*PHI20;

% levitation forces
L1 = (mu0*Nr^2*Al/4)*(ir1/zFL)^2;
L2 = (mu0*Nr^2*Al/4)*(ir2/zCL)^2;
L3 = (mu0*Nr^2*Al/4)*(ir3/zRL)^2;
L4 = (mu0*Nr^2*Al/4)*(ir4/zFR)^2;
L5 = (mu0*Nr^2*Al/4)*(ir5/zCR)^2;
L6 = (mu0*Nr^2*Al/4)*(ir6/zRR)^2;

%distribuzioni delle forze sul modello elastico
lf=(L1+L2-((rho*Ss*Cl/2)*(vx)^2+m*radius*(thetadot^2)*sin(phi)+m*g*cos(phi))/3)/l;
lc=(L3+L4-((rho*Ss*Cl/2)*(vx)^2+m*radius*(thetadot^2)*sin(phi)+m*g*cos(phi))/3)/l;
lr=(L5+L6-((rho*Ss*Cl/2)*(vx)^2+m*radius*(thetadot^2)*sin(phi)+m*g*cos(phi))/3)/l;

% diff equations
z_dot    = zdot;
%togliere cos(theta) sotto. sostituisci m*radius
%***********************************************************************
zdot_dot = (L1 + L2 + L3 + L4 + L5 + L6 - (rho*Ss*Cl/2)*(xdot_ref)^2 - m*g*cos(phi)*cos(theta) - m*radius*(thetadot^2)*sin(phi))/m;% + xdot_ref*thetadot - vy*wphi

phi_dot    = phidot;
phidot_dot = (yB*(L1 + L2 + L3 - L4 - L5 - L6) + zB*(G1 - G4 + G2 - G5 + G3 - G6))/Jx;

theta_dot    = thetadot;
thetadot_dot = (xB*(-L1 - L4 + L3 + L6) -zA*(rho*Ss*Cd/2)*(xdot_ref)^2 + xA*(rho*Ss*Cl/2)*(xdot_ref)^2 +zB*(P1+P2+P3+P4+P5+P6))/Jy;

ir1_dot = (ir1/zFL)*zFL_dot -(2/(mu0*(Nr^2)*Al))*zFL*(Rz*ir1-Vr1);
ir2_dot = (ir2/zCL)*zCL_dot -(2/(mu0*(Nr^2)*Al))*zCL*(Rz*ir2-Vr2);
ir3_dot = (ir3/zRL)*zRL_dot -(2/(mu0*(Nr^2)*Al))*zRL*(Rz*ir3-Vr3);
ir4_dot = (ir4/zFR)*zFR_dot -(2/(mu0*(Nr^2)*Al))*zFR*(Rz*ir4-Vr4);
ir5_dot = (ir5/zCR)*zCR_dot -(2/(mu0*(Nr^2)*Al))*zCR*(Rz*ir5-Vr5);
ir6_dot = (ir6/zRR)*zRR_dot -(2/(mu0*(Nr^2)*Al))*zRR*(Rz*ir6-Vr6);


%EQUAZIONI DEL MODELLO ELASTICO (per ora ne ho scritte solo due, poi vediamo)
qv_dot=vqv;
vqv_dot=-qv*w^2+lf*PHIf+lc*PHIc+lr*PHIr;

%differential equation vector
ode = [z_dot zdot_dot phi_dot phidot_dot theta_dot thetadot_dot ir1_dot ir2_dot ir3_dot ir4_dot ir5_dot ir6_dot qv_dot vqv_dot]';

% Linearization

% references
z_ref = 1.04;
zdot_ref = 0;
phi_ref = 0;
phidot_ref = 0;
theta_ref = 0;
thetadot_ref = 0;
qv_ref = 0;
vqv_ref = 0;
zi_ref = h - z_ref +zi_rip;


iri_ref = sqrt( (1/(6*(mu0*Nr^2*Al/4)))*((rho*Ss*Cl/2)*(xdot_ref)^2 + m*g*cos(phi_ref) + m*radius*(thetadot_ref^2)*sin(phi_ref))*zi_ref^2 -m*xdot_ref*thetadot_ref);
L2_ref = (mu0*Nr^2*Al/4)*(iri_ref/zi_ref)^2;
L5_ref = L2_ref;
L1_ref = -(1/2)*L2_ref -(1/4)*((zA/xB)*(rho*Ss*Cd/2)*xdot_ref^2 - (xA/xB)*(rho*Ss*Cl/2)*xdot_ref^2 - (zB/xB)*(P1+P2+P3+P4+P5+P6)) + (1/4)*(rho*Ss*Cl/2)*xdot_ref^2 + (1/4)*m*g*cos(phi_ref)*cos(theta_ref); %-(1/4)*xdot_ref*thetadot_ref;
L4_ref = L1_ref;
L3_ref = (1/2)*((zA/xB)*(rho*Ss*Cd/2)*xdot_ref^2 - (xA/xB)*(rho*Ss*Cl/2)*xdot_ref^2 - (zB/xB)*(P1+P2+P3+P4+P5+P6)+2*L1_ref);
L6_ref = L3_ref;
ir1_ref = sqrt( 1/(mu0*Nr^2*Al/4)*(L1_ref)*zi_ref^2 );
ir2_ref = sqrt( 1/(mu0*Nr^2*Al/4)*(L2_ref)*zi_ref^2 );
ir3_ref = sqrt( 1/(mu0*Nr^2*Al/4)*(L3_ref)*zi_ref^2 );
ir4_ref = sqrt( 1/(mu0*Nr^2*Al/4)*(L4_ref)*zi_ref^2 );
ir5_ref = sqrt( 1/(mu0*Nr^2*Al/4)*(L5_ref)*zi_ref^2 );
ir6_ref = sqrt( 1/(mu0*Nr^2*Al/4)*(L6_ref)*zi_ref^2 );


Vr1_ref = Rz*ir1_ref;
Vr2_ref = Rz*ir2_ref;
Vr3_ref = Rz*ir3_ref;
Vr4_ref = Rz*ir4_ref;
Vr5_ref = Rz*ir5_ref;
Vr6_ref = Rz*ir6_ref;

Vref_z=[Vr1_ref; Vr2_ref; Vr3_ref; Vr4_ref; Vr5_ref; Vr6_ref];

Az = jacobian(ode,states);
Az = subs(Az, {z zdot phi phidot theta thetadot ir1 ir2 ir3 ir4 ir5 ir6 qv vqv Vr1 Vr2 Vr3 Vr4 Vr5 Vr6}, {z_ref zdot_ref phi_ref phidot_ref theta_ref thetadot_ref ir1_ref ir2_ref ir3_ref ir4_ref ir5_ref ir6_ref qv_ref vqv_ref Vr1_ref Vr2_ref Vr3_ref Vr4_ref Vr5_ref Vr6_ref});
Az = eval(Az);
%A = vpa(A,5)

Bz = jacobian(ode,inputs);
Bz = subs(Bz,{theta phi z qv vqv},{theta_ref phi_ref z_ref qv_ref vqv_ref});
Bz = eval(Bz);
%B = vpa(B, 5)

Cz = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 ;
    0 0 1 0 0 0 0 0 0 0 0 0 0 0 ;
    0 0 0 0 1 0 0 0 0 0 0 0 0 0 ;
    0 0 0 0 0 0 1 0 0 0 0 0 0 0 ;
    0 0 0 0 0 0 0 1 0 0 0 0 0 0 ;
    0 0 0 0 0 0 0 0 1 0 0 0 0 0 ;
    0 0 0 0 0 0 0 0 0 1 0 0 0 0 ;
    0 0 0 0 0 0 0 0 0 0 1 0 0 0 ;
    0 0 0 0 0 0 0 0 0 0 0 1 0 0 ];

Dz = zeros(9,6);

contrz = ctrb(Az,Bz);
rz = rank(contrz);


% LQR control
i_cost=10000000000;
% Q = diag([1 1 1 1 1 1 i_cost i_cost i_cost i_cost i_cost i_cost]);
Qz = diag([1 0.1 1 0.1 1 1 i_cost i_cost i_cost i_cost i_cost i_cost 1 1]);
Rrz = 1000000000*eye(6);
Rrkz= 100000000*eye(9);
Kz = lqr(Az,Bz,Qz,Rrz);

Zref = [z_ref zdot_ref phi_ref phidot_ref theta_ref thetadot_ref ir1_ref ir2_ref ir3_ref ir4_ref ir5_ref ir6_ref qv_ref vqv_ref].';

%initial conditions for nonlinear sim
phi_init=0;
theta_init=0;
qv_init = 0;
vqv_init = 0;

Z_init = [h; 0; phi_init; 0; theta_init; 0; 0;0;0;0;0;0;0;0];