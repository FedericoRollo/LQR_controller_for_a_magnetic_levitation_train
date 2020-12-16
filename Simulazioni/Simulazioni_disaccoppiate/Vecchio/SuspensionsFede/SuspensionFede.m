clear all; clc;

%Constants
xB=1;
yB=1;
PHI=1;
PHIf=1;
PHIc=1;
PHIr=1;
mu0=1;
Nr=1;
Al=1;
rho=1;
S=1;
Cd=1;
Cl=1;
g=9.81;
R=1;
m=1;
zB=1;
xA=1;
zA=1;
Itheta=1;
alpha=1;
P1=1; 
P2=1;
P3=1; 
P4=1; 
P5=1; 
P6=1;
G1=1;
G2=1;
G3=1;
G4=1;
G5=1;
G6=1;
vx=1;
Iphi=1;
w=1;
l=1;

%references
z_ref       = 1;
vz_ref      = 0;
phi_ref     = 0;
wphi_ref    = 0;
theta_ref   = 0;
wtheta_ref  = 0;
ir1_ref     = rand(1);
ir2_ref     = rand(1);
ir3_ref     = rand(1);
ir4_ref     = rand(1);
ir5_ref     = rand(1);
ir6_ref     = rand(1);
qv_ref      = rand(1);
vqv_ref     = rand(1);

%equilibrum point
xe = transpose([z_ref vz_ref phi_ref wphi_ref theta_ref wtheta_ref ir1_ref ir2_ref ir3_ref ir4_ref ir5_ref ir6_ref qv_ref vqv_ref]);

% xe = transpose([z_ref vz_ref phi_ref wphi_ref theta_ref wtheta_ref ir1_ref ir2_ref ir3_ref ir4_ref ir5_ref ir6_ref]);

Vr1_ref = rand(1);
Vr2_ref = rand(1);
Vr3_ref = rand(1);
Vr4_ref = rand(1);
Vr5_ref = rand(1);
Vr6_ref = rand(1);

%% Differential equation definition
syms z vz phi wphi theta wtheta ir1 ir2 ir3 ir4 ir5 ir6 qv vqv Vr1 Vr2 Vr3 Vr4 Vr5 Vr6 

% syms z vz phi wphi theta wtheta ir1 ir2 ir3 ir4 ir5 ir6 Vr1 Vr2 Vr3 Vr4 Vr5 Vr6 

stati = transpose([z vz phi wphi theta wtheta ir1 ir2 ir3 ir4 ir5 ir6 qv vqv]);
% stati = transpose([z vz phi wphi theta wtheta ir1 ir2 ir3 ir4 ir5 ir6]);
input = transpose([Vr1 Vr2 Vr3 Vr4 Vr5 Vr6]);


%spostamenti verticali dei rotori e dei sensori
deltaz1=-z-phi*yB+theta*xB-qv*PHI;
deltaz2=-z+phi*yB+theta*xB-qv*PHI;
deltaz3=-z-phi*yB-qv*PHI;
deltaz4=-z+phi*yB-qv*PHI;
deltaz5=-z-phi*yB-theta*xB-qv*PHI;
deltaz6=-z+phi*yB-theta*xB-qv*PHI;

deltazdot1=-vz-wphi*yB+wtheta*xB-vqv*PHI;
deltazdot2=-vz+wphi*yB+wtheta*xB-vqv*PHI;
deltazdot3=-vz-wphi*yB-vqv*PHI;
deltazdot4=-vz+wphi*yB-vqv*PHI;
deltazdot5=-vz-wphi*yB-wtheta*xB-vqv*PHI;
deltazdot6=-vz+wphi*yB-wtheta*xB-vqv*PHI;

%forze di levitazione
L1=(mu0*Nr^2*Al/4)*(ir1/deltaz1)^2;
L2=(mu0*Nr^2*Al/4)*(ir2/deltaz2)^2;
L3=(mu0*Nr^2*Al/4)*(ir3/deltaz3)^2;
L4=(mu0*Nr^2*Al/4)*(ir4/deltaz4)^2;
L5=(mu0*Nr^2*Al/4)*(ir5/deltaz5)^2;
L6=(mu0*Nr^2*Al/4)*(ir6/deltaz6)^2;

%distribuzioni delle forze sul modello elastico
lf=(L1+L2-((rho*S*Cl/2)*(vx)^2+m*R*(wtheta^2)*sin(phi)+m*g*cos(phi))/3)/l;
lc=(L3+L4-((rho*S*Cl/2)*(vx)^2+m*R*(wtheta^2)*sin(phi)+m*g*cos(phi))/3)/l;
lr=(L5+L6-((rho*S*Cl/2)*(vx)^2+m*R*(wtheta^2)*sin(phi)+m*g*cos(phi))/3)/l;

%z axis movement
zdot=vz;
vzdot=(L1+L2+L3+L4+L5+L6-(rho*S*Cl/2)*(vz)^2-m*R*(wtheta^2)*sin(phi)-m*g*cos(phi))/m;
%roll angle (around x)
phidot=wphi;
wphidot=(yB*(L1-L2+L3-L4+L5-L6)+zB*(G1-G2+G3-G4+G5-G6))/Iphi;
%pitch angle (around y)
thetadot=wtheta;
wthetadot=(xB*(-L1-L2+L5+L6)+zB*(P1+P2+P3+P4+P5+P6)-zA*(rho*S*Cd/2)*(vx)^2+xA*(rho*S*Cl/2)*(vx)^2)/Itheta;

%LKT DEI ROTORI
irdot1=ir1*deltazdot1/deltaz1-(2/(mu0*Nr*alpha))*deltaz1*(R*ir1-Vr1);
irdot2=ir2*deltazdot2/deltaz2-(2/(mu0*Nr*alpha))*deltaz2*(R*ir2-Vr2);
irdot3=ir3*deltazdot3/deltaz3-(2/(mu0*Nr*alpha))*deltaz3*(R*ir3-Vr3);
irdot4=ir4*deltazdot4/deltaz4-(2/(mu0*Nr*alpha))*deltaz4*(R*ir4-Vr4);
irdot5=ir5*deltazdot5/deltaz5-(2/(mu0*Nr*alpha))*deltaz5*(R*ir5-Vr5);
irdot6=ir6*deltazdot6/deltaz6-(2/(mu0*Nr*alpha))*deltaz6*(R*ir6-Vr6);

%EQUAZIONI DEL MODELLO ELASTICO (per ora ne ho scritte solo due, poi vediamo)
qvdot=vqv;
vqvdot=-qv*w^2+lf*PHIf+lc*PHIc+lr*PHIr;

%differential equation vector
ode = transpose([zdot vzdot phidot wphidot thetadot wthetadot irdot1 irdot2 irdot3 irdot4 irdot5 irdot6 qvdot vqvdot]);

% ode = transpose([zdot vzdot phidot wphidot thetadot wthetadot irdot1 irdot2 irdot3 irdot4 irdot5 irdot6]);

%% Linearization
a = jacobian(ode,stati);
b = jacobian(ode,input);

z       = z_ref;
vz      = vz_ref;
phi     = phi_ref;
wphi    = wphi_ref;
theta   = theta_ref;
wtheta  = wtheta_ref;
ir1     = ir1_ref;
ir2     = ir2_ref;
ir3     = ir3_ref;
ir4     = ir4_ref;
ir5     = ir5_ref;
ir6     = ir6_ref;
qv      = qv_ref;
vqv     = vqv_ref;

Vr1     = Vr1_ref;
Vr2     = Vr2_ref;
Vr3     = Vr3_ref;
Vr4     = Vr4_ref;
Vr5     = Vr5_ref;
Vr6     = Vr6_ref;

A = eval(a);
B = eval(b);

contr = ctrb(A,B);
r = rank(contr);


%% LQR control

Q = diag([100 1 10 1 10 1 1 1 1 1 1 1 1 1 ]);
% Q = diag([100 1 10 1 10 1 1 1 1 1 1 1 ]);
R = 1*eye(6);

K = lqr(A,B,Q,R);

