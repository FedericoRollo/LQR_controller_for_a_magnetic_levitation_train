%% HEADWAY motion
% we are considering 2 stators and 6 rotors, the 6 rotors have a DC
% excitation considered for the forward motion in the nominal value. In
% this way we decouple its motion from the others
clear all;
clc;

%% constants
%********************LE COSTANTI VANNO CABIATEEEEE*************************
% una volta cambiate le costanti vanno risistemati anche i vari valori per
% il controllo ottimo perchè non recupera automaticamente le incertezze!!

%Parametri brutti

% z=1; % nominal z
% i_e=1; % nominal current which provide nominal z
% R=1;
% Ld=1;
% Lq=1;
% M=10;%mutual inductance stator rotor
% tau=1;
% p=3;% poli!!!!!!!!!!!!!!!!!!!!!!
% m=1;%mass 
% rho=1;
% S=1;
% Cd=1;

%Parametri Adelina

z=0.05; % nominal z
i_e=100; % nominal current which provides nominal z
R=1.04; %resistance of the stator
Ld=0.06;
Lq=0.06;
M=0.1;%mutual inductance stator rotor
tau=1; %pole pitch
p=3;% poles
m=2340;%mass 
rho=1.225; %aerodynamic coefficient
S=10; %aerodynamic coefficient
Cd=0.8; %aerodynamic coefficient

%references THEY'RE USED for linearization and ALSO IN THE SIMULATION SO THEY ARE IMPORTANT
iq1_ref = 0;
id1_ref = 0;
iq2_ref = 0;
id2_ref = 0;
v_ref   = 100;

xe = [iq1_ref; id1_ref; iq2_ref; id2_ref; v_ref]; %needed only for the
% linear simulation while non linear simulation need a speed profile
% otherwise it diverges

%IMPORTANT FOR SIMULATION: also a profile for iq can be used in order to
%obtain a smooth behaviour for the current and to attenuate in some way the
%input voltages

Vq1_ref=rand(1);
Vd1_ref=rand(1);
Vq2_ref=rand(1);
Vd2_ref=rand(1);



%% stators currents

syms iq1 iq2 id1 id2 v Vq1 Vq2 Vd1 Vd2

states = [iq1; id1; iq2; id2; v];
inputs = [Vq1; Vd1; Vq2; Vd2];

iq1_dot = -(R/Lq)*iq1 + (Ld/Lq)*((M/Ld)*i_e + id1)*v*pi/tau + Vq1/Lq ;
iq2_dot = -(R/Lq)*iq2 + (Ld/Lq)*((M/Ld)*i_e + id2)*v*pi/tau + Vq2/Lq ;

id1_dot = -(R/Ld)*id1 - (Lq/Ld)*(iq1)*v*pi/tau + Vd1/Ld ;
id2_dot = -(R/Ld)*id2 - (Lq/Ld)*(iq2)*v*pi/tau + Vd2/Ld ;

Fx1 = (3/2)*(pi/tau)*(Ld*((M/Ld)*i_e) + (Ld-Lq)*id1)*iq1;
Fx2 = (3/2)*(pi/tau)*(Ld*((M/Ld)*i_e) + (Ld-Lq)*id2)*iq2;

v_dot = (1/m)*(p*Fx1 + p*Fx2 -(rho*S*Cd/2)*v^2);

%differential equation vector
ode = [iq1_dot; id1_dot; iq2_dot; id2_dot; v_dot];

%% Linearization

A= jacobian(ode, states);
B= jacobian(ode, inputs);

%for LQR C and D are useless because we know all the state but they can be
%used for LQG if C changes
C=eye(5);
D=[zeros(4); 0 0 0 0];

iq1 = iq1_ref;
id1 = id1_ref;
iq2 = iq2_ref;
id2 = id2_ref;
v   = v_ref;

Vq1 = Vq1_ref; 
Vd1 = Vd1_ref;
Vq2 = Vq2_ref;
Vd2 = Vd2_ref;

Fx1_ref = ((rho*S*Cd/2)*v_ref^2)/(2*p);
Fx2_ref = Fx1_ref;

iq1_ref = (Fx1_ref)/((3/2)*(pi/tau)*(Ld*(M/Ld)*i_e));
iq2_ref = (Fx2_ref)/((3/2)*(pi/tau)*(Ld*(M/Ld)*i_e));

Vq1_ref = R*iq1_ref-Ld*((M/Ld)*i_e+id1_ref)*v_ref*(pi/tau);
Vq2_ref = R*iq2_ref-Ld*((M/Ld)*i_e+id2_ref)*v_ref*(pi/tau);

Vd1_ref = R*id1_ref+Lq*(iq1_ref)*v_ref*pi/tau ;
Vd2_ref = R*id2_ref+Lq*(iq2_ref)*v_ref*pi/tau ;

A=eval(A);
B=eval(B);

%the values must be choosen also depending on previous values!!
Q=[10 0 0 0 0; 0 10 0 0 0; 0 0 10 0 0; 0 0 0 10 0; 0 0 0 0 100];
R=1*eye(4);

xe = [iq1_ref; id1_ref; iq2_ref; id2_ref; v_ref];  

K=lqr(A,B,Q,R);

%Linear simulation
open('HeadwayLinearFede.slx');
sim('HeadwayLinearFede.slx');

%Non-Linear Simulation
open('HeadwayNonLinearFede.slx');
sim('HeadwayNonLinearFede.slx');


