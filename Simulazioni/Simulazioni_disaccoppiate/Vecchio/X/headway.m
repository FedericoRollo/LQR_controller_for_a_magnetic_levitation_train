%% HEADWAY motion
% we are considering 2 stators and 6 rotors, the 6 rotors have a DC
% excitation considered for the forward motion in the nominal value. In
% this way we decouple its motion from the others

clear all
clc

%% constants
z=0.05; % nominal z

%motor constants
i_e=100; % nominal current which provides nominal z (it MUST be replaced with the REAL NOMINAL VALUE)
R=1.04; %resistance of the stator
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

iq1_ref=0;
id1_ref=0;
iq2_ref=0;
id2_ref=0;

xdot_ref=430/3.6; %m/s

Vq1_ref=1;
Vd1_ref=1;
Vq2_ref=1;
Vd2_ref=1;



%% stators currents

syms iq1 iq2 id1 id2 xdot Vq1 Vq2 Vd1 Vd2

states = [iq1; id1; iq2; id2; xdot];
inputs = [Vq1; Vd1; Vq2; Vd2];

iq1_dot = -(R/Lq)*iq1 + (Ld/Lq)*((M/Ld)*i_e + id1)*xdot*pi/tau + Vq1/Lq ;
iq2_dot = -(R/Lq)*iq2 + (Ld/Lq)*((M/Ld)*i_e + id2)*xdot*pi/tau + Vq2/Lq ;

id1_dot = -(R/Ld)*id1 - (Lq/Ld)*(iq1)*xdot*pi/tau + Vd1/Ld ;
id2_dot = -(R/Ld)*id2 - (Lq/Ld)*(iq2)*xdot*pi/tau + Vd2/Ld ;

Fx1 = (3/2)*(pi/tau)*(Ld*((M/Ld)*i_e) + (Ld-Lq)*id1)*iq1;
Fx2 = (3/2)*(pi/tau)*(Ld*((M/Ld)*i_e) + (Ld-Lq)*id2)*iq2;

xdot_dot = (1/m)*(p*Fx1 + p*Fx2 - (rho*Sa*Cd/2)*xdot^2);

ode = [iq1_dot; id1_dot; iq2_dot; id2_dot; xdot_dot];

A=jacobian(ode, states);

%references for currents, related to speed
Fx1_ref = ((rho*Sa*Cd/2)*xdot_ref^2)/(2*p);
Fx2_ref = Fx1_ref;

iq1_ref = (Fx1_ref)/((3/2)*(pi/tau)*(Ld*(M/Ld)*i_e));
iq2_ref = (Fx2_ref)/((3/2)*(pi/tau)*(Ld*(M/Ld)*i_e));

Vq1_ref = R*iq1_ref-Ld*((M/Ld)*i_e+id1_ref)*xdot_ref*(pi/tau);
Vq2_ref = R*iq2_ref-Ld*((M/Ld)*i_e+id2_ref)*xdot_ref*(pi/tau);

Vd1_ref = R*id1_ref+Lq*(iq1_ref)*xdot_ref*pi/tau ;
Vd2_ref = R*id2_ref+Lq*(iq2_ref)*xdot_ref*pi/tau ;


A=subs(A, {iq1 id1 iq2 id2 xdot}, {iq1_ref id1_ref iq2_ref id2_ref xdot_ref});
A=eval(A);
B= jacobian(ode, inputs);

B=subs(B, {Vq1 Vd1 Vq2 Vd2}, {Vq1_ref Vd1_ref Vq2_ref Vd2_ref});
B=eval(B);
Q=[1000 0 0 0 0; 0 100000 0 0 0; 0 0 1000 0 0; 0 0 0 100000 0; 0 0 0 0 10];
Rr=[10 0 0 0; 0 10 0 0; 0 0 10 0; 0 0 0 10];

x_ref = [iq1_ref; id1_ref; iq2_ref; id2_ref; xdot_ref];
u_ref = [Vq1_ref; Vd1_ref; Vq2_ref; Vd2_ref];

[K,S,e]=lqr(A,B,Q,Rr);