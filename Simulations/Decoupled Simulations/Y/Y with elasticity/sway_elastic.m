%% SWAY MOTION

clear all
clc

% Constants


vx=430/3.6;
xB=8;
yB=2;
mu0=1.25*10^-6;
Ng=500;
Ag=1;
m=60000;
g=9.8;
Ra=0;
R=100;
Jz=2853*10^3;
phi=0;
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



%% ------------------------STATES STATEMENT--------------------------------
syms y vy psi wpsi ql vql ig1 ig2 ig3 ig4 ig5 ig6 

states = [ y vy psi wpsi ig1 ig2 ig3 ig4 ig5 ig6 ql vql ];
% states = [ y vy psi wpsi ig1 ig2 ig3 ig4 ig5 ig6 ];

%% -----------------------INPUT STATEMENT----------------------------------
syms Vg1 Vg2 Vg3 Vg4 Vg5 Vg6 
input = [Vg1 Vg2 Vg3 Vg4 Vg5 Vg6];

%% Differential equations
%1;2;3 -> right electromagnets, 4;5;6 -> left
%2;4;6                          1;3;5
% horizontal displacement

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
gf=(G4-G1 + (m*g*sin(phi)-m*Ra*cos(phi)*wpsi^2)/3)/l;
gc=(G5-G2 +(m*g*sin(phi)-m*Ra*cos(phi)*wpsi^2)/3)/l;
gr=(G6-G3 +(m*g*sin(phi)-m*Ra*cos(phi)*wpsi^2)/3)/l;


%----------differential equations-----------

%z - I'm using only this angle for two reason: because the other two are
%considered in the suspentions and also because this angles could not
%affect the distance from the guidance elector magnet if they're at the
%same height of the center of mass
%1;2;3 -> right electromagnets, 4;5;6 -> left
%2;4;6                          1;3;5
psidot  = wpsi;
wpsidot = ((xB*(-G4+G1+G6-G3)+yB*(-P1+P2-P3+P4-P5+P6))/Jz); %(Jy-Jz)*tethadot*psidot -> si toglie
ydot  = vy;
vydot=((-G4+G1-G5+G2-G6+G3+m*g*sin(phi)-m*Ra*(wpsi^2)*cos(phi))/m) - vx*psidot;% + vz*phidot -> si toglie


%equazioni elastiche
qldot=vql;
vqldot=-ql*w^2+gf*PHIf+gc*PHIc+gr*PHIr;


%LKT DEGLI ELETTROAMGNETI DELLA GUIDANCE
igdot1=ig1*deltaydot1/deltay1-(2/(mu0*(Ng^2)*Ag))*deltay1*(R*ig1-Vg1);
igdot2=ig2*deltaydot2/deltay2-(2/(mu0*(Ng^2)*Ag))*deltay2*(R*ig2-Vg2);
igdot3=ig3*deltaydot3/deltay3-(2/(mu0*(Ng^2)*Ag))*deltay3*(R*ig3-Vg3);
igdot4=ig4*deltaydot4/deltay4-(2/(mu0*(Ng^2)*Ag))*deltay4*(R*ig4-Vg4);
igdot5=ig5*deltaydot5/deltay5-(2/(mu0*(Ng^2)*Ag))*deltay5*(R*ig5-Vg5);
igdot6=ig6*deltaydot6/deltay6-(2/(mu0*(Ng^2)*Ag))*deltay6*(R*ig6-Vg6);


%diffetential equation
ode = [ydot; vydot; psidot; wpsidot; igdot1; igdot2; igdot3; igdot4; igdot5; igdot6; qldot; vqldot];
% ode = [ydot; vydot; psidot; wpsidot; igdot1; igdot2; igdot3; igdot4; igdot5; igdot6];

%% Linearization

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
Vg_ref = R*ig_ref;
Vg1 = Vg_ref;
Vg2 = Vg_ref;
Vg3 = Vg_ref;
Vg4 = Vg_ref;
Vg5 = Vg_ref;
Vg6 = Vg_ref;

Vref = [ Vg1; Vg2; Vg3; Vg4; Vg5; Vg6];
%compute A B matricies
A=eval(a);
% A=subs(a,[y vy psi wpsi ig1 ig2 ig3 ig4 ig5 ig6 ql vql Vg1 Vg2 Vg3 Vg4 Vg5 Vg6], [y_ref vy_ref psi_ref wpsi_ref ig_ref ig_ref ig_ref ig_ref ig_ref ig_ref ql_ref vql_ref Vg_ref Vg_ref Vg_ref Vg_ref Vg_ref Vg_ref]);
B=eval(b);
% B=subs(b,[y vy psi wpsi ig1 ig2 ig3 ig4 ig5 ig6 ql vql Vg1 Vg2 Vg3 Vg4 Vg5 Vg6], [y_ref vy_ref psi_ref wpsi_ref ig_ref ig_ref ig_ref ig_ref ig_ref ig_ref ql_ref vql_ref Vg_ref Vg_ref Vg_ref Vg_ref Vg_ref Vg_ref]);

Co = ctrb(A,B); % controllability matrix
r = rank(Co);

%{
% set Q and R matricies values respectively for states and input weights
weight_y      = 1;
weight_vy     = 0.0001;
weight_psi    = 0.0001;
weight_wpsi   = 0.0001;
weight_ql     = 0.0001;
weight_vql    = 0.0001;
weight_ig1    = 0.0001;
weight_ig2    = 0.0001;
weight_ig3    = 0.0001;
weight_ig4    = 0.0001;
weight_ig5    = 0.0001;
weight_ig6    = 0.0001;

% Q = diag([ weight_y weight_vy weight_psi weight_wpsi weight_ql weight_vql weight_ig1 weight_ig2 weight_ig3 weight_ig4 weight_ig5 weight_ig6]);
Q = diag([ weight_y weight_vy weight_psi weight_wpsi weight_ig1 weight_ig2 weight_ig3 weight_ig4 weight_ig5 weight_ig6]);

weight_Vg1 = 0.0001;
weight_Vg2 = 0.0001;
weight_Vg3 = 0.0001;
weight_Vg4 = 0.0001;
weight_Vg5 = 0.0001;
weight_Vg6 = 0.0001;

R = diag([weight_Vg1 weight_Vg2 weight_Vg3 weight_Vg4 weight_Vg5 weight_Vg6]);
%}

% Q = diag([100 100 10 1 10 10 10 10 10 10]);
Q = diag([1000 0.1 10 10 100000 100000 100000 100000 100000 100000 0.1 0.1]);
Rr = 1*eye(6);

K = lqr(A,B,Q,Rr);

%INITIAL VALUES

%set states inital values
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


x_ref = [y_ref; vy_ref; psi_ref; wpsi_ref; ig_ref; ig_ref; ig_ref; ig_ref; ig_ref; ig_ref; ql_ref; vql_ref];
