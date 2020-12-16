
%% Pump motion
% vertical motion (along z axes) is here analyzed. current 1 -> right front; 2 -> right center; 3-> right rear; 4 -> left front; 5 -> left center; 6 -> left rear   
clear all
clc
%% Constants
h = 1; %constant, height of z at rest
zi_rip = 0.1; %zi constant of zi when the train is at rest: graphically z(h)|zi(-zi_rest), where |->axes
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
R=100;
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
vx=100;
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

%% Differential equation
 
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
lf=(L1+L2-((rho*Ss*Cl/2)*(vx)^2+m*R*(thetadot^2)*sin(phi)+m*g*cos(phi))/3)/l;
lc=(L3+L4-((rho*Ss*Cl/2)*(vx)^2+m*R*(thetadot^2)*sin(phi)+m*g*cos(phi))/3)/l;
lr=(L5+L6-((rho*Ss*Cl/2)*(vx)^2+m*R*(thetadot^2)*sin(phi)+m*g*cos(phi))/3)/l;


% diff equations
z_dot    = zdot;
zdot_dot = (L1 + L2 + L3 + L4 + L5 + L6 - (rho*Ss*Cl/2)*(xdot_ref)^2 - m*g*cos(phi)*cos(theta) - m*radius*(thetadot^2)*sin(phi))/m;% +xdot_ref*thetadot

phi_dot    = phidot;
phidot_dot = (yB*(L1 + L2 + L3 - L4 - L5 - L6) + zB*(G1 - G4 + G2 - G5 + G3 - G6))/Jx;

theta_dot    = thetadot;
thetadot_dot = (xB*(-L1 - L4 + L3 + L6) -zA*(rho*Ss*Cd/2)*(xdot_ref)^2 + xA*(rho*Ss*Cl/2)*(xdot_ref)^2 +zB*(P1+P2+P3+P4+P5+P6))/Jy;

ir1_dot = (ir1/zFL)*zFL_dot -(2/(mu0*(Nr^2)*Al))*zFL*(R*ir1-Vr1);
ir2_dot = (ir2/zCL)*zCL_dot -(2/(mu0*(Nr^2)*Al))*zCL*(R*ir2-Vr2);
ir3_dot = (ir3/zRL)*zRL_dot -(2/(mu0*(Nr^2)*Al))*zRL*(R*ir3-Vr3);
ir4_dot = (ir4/zFR)*zFR_dot -(2/(mu0*(Nr^2)*Al))*zFR*(R*ir4-Vr4);
ir5_dot = (ir5/zCR)*zCR_dot -(2/(mu0*(Nr^2)*Al))*zCR*(R*ir5-Vr5);
ir6_dot = (ir6/zRR)*zRR_dot -(2/(mu0*(Nr^2)*Al))*zRR*(R*ir6-Vr6);


%EQUAZIONI DEL MODELLO ELASTICO (per ora ne ho scritte solo due, poi vediamo)
qv_dot=vqv;
vqv_dot=-qv*w^2+lf*PHIf+lc*PHIc+lr*PHIr;


%differential equation vector
ode = [z_dot zdot_dot phi_dot phidot_dot theta_dot thetadot_dot ir1_dot ir2_dot ir3_dot ir4_dot ir5_dot ir6_dot qv_dot vqv_dot]';
% ode = [z_dot zdot_dot phi_dot phidot_dot theta_dot thetadot_dot ir1_dot ir2_dot ir3_dot ir4_dot ir5_dot ir6_dot]';
% ode = transpose([zdot vzdot phidot wphidot thetadot wthetadot irdot1 irdot2 irdot3 irdot4 irdot5 irdot6]);

%% Linearization

% references
z_ref = 1.09;
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


Vr1_ref = R*ir1_ref;
Vr2_ref = R*ir2_ref;
Vr3_ref = R*ir3_ref;
Vr4_ref = R*ir4_ref;
Vr5_ref = R*ir5_ref;
Vr6_ref = R*ir6_ref;

Vref=[Vr1_ref; Vr2_ref; Vr3_ref; Vr4_ref; Vr5_ref; Vr6_ref;]

A = jacobian(ode,states);
A = subs(A, {z zdot phi phidot theta thetadot ir1 ir2 ir3 ir4 ir5 ir6 qv vqv Vr1 Vr2 Vr3 Vr4 Vr5 Vr6}, {z_ref zdot_ref phi_ref phidot_ref theta_ref thetadot_ref ir1_ref ir2_ref ir3_ref ir4_ref ir5_ref ir6_ref qv_ref vqv_ref Vr1_ref Vr2_ref Vr3_ref Vr4_ref Vr5_ref Vr6_ref});
A = eval(A)
%A = vpa(A,5)

B = jacobian(ode,inputs);
B = subs(B,{theta phi z qv vqv},{theta_ref phi_ref z_ref qv_ref vqv_ref});
B = eval(B) 
%B = vpa(B, 5)

contr = ctrb(A,B);
r = rank(contr)


%% LQR control
i_cost=100000;
% Q = diag([1 1 1 1 1 1 i_cost i_cost i_cost i_cost i_cost i_cost]);
Q = diag([1 0.1 1 0.1 1 1 i_cost i_cost i_cost i_cost i_cost i_cost 1 1]);
Rr = 1000*eye(6);

K = lqr(A,B,Q,Rr);

xe = [z_ref zdot_ref phi_ref phidot_ref theta_ref thetadot_ref ir1_ref ir2_ref ir3_ref ir4_ref ir5_ref ir6_ref qv_ref vqv_ref].';
%xe = eval(xe)
%initial conditions for nonlinear sim
phi_init=0;
theta_init=0;
qv_init = 0;
vqv_init = 0;