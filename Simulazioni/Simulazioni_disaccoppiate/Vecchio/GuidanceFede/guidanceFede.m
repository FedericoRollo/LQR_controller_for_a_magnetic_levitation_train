clear all; clc;
%------------------------STATES STATEMENT--------------------------------

syms y vy psi wpsi ql vql ig1 ig2 ig3 ig4 ig5 ig6 

% stati = [ y vy psi wpsi ql vql ig1 ig2 ig3 ig4 ig5 ig6 ];
stati = [ y vy psi wpsi ig1 ig2 ig3 ig4 ig5 ig6 ];

%-----------------------INPUT STATEMENT----------------------------------
syms Vg1 Vg2 Vg3 Vg4 Vg5 Vg6 
input = [Vg1 Vg2 Vg3 Vg4 Vg5 Vg6];

%Constants
xB=1;
yB=1;
mu0=1;
Ng=1;
Ag=1;
m=1;
g=9.8;
Ra=1;
Ipsi=1;
w=1;
alpha=1;
phi=1;
wtheta=1;
P1=1;
P2=1;
P3=1;
P4=1;
P5=1;
P6=1;



%spostamenti orizzontali degli elettromagneti e dei sensori
deltay2=-y-psi*xB%-ql*PHI;
deltay4=-y%-ql*PHI;
deltay6=-y+psi*xB%-ql*PHI;
deltay1=-deltay2;
deltay3=-deltay4;
deltay5=-deltay6;

deltaydot2=-vy-wpsi*xB%-vql*PHI;
deltaydot4=-vy%-vql*PHI;
deltaydot6=-vy+wpsi*xB%-vql*PHI;
deltaydot1=-deltaydot2;
deltaydot3=-deltaydot4;
deltaydot5=-deltaydot6;

% %distribuzioni delle forze sul modello elastico
% gf=(G1-G2+(m*g*sin(phi)-m*Ra*cos(phi)*wtheta^2)/3)/l;
% gc=(G3-G4+(m*g*sin(phi)-m*Ra*cos(phi)*wtheta^2)/3)/l;
% gr=(G5-G6+(m*g*sin(phi)-m*Ra*cos(phi)*wtheta^2)/3)/l;

%forse della guidance
G1=(mu0*Ng^2*Ag/4)*(ig1/deltay1)^2;
G2=(mu0*Ng^2*Ag/4)*(ig2/deltay2)^2;
G3=(mu0*Ng^2*Ag/4)*(ig3/deltay3)^2;
G4=(mu0*Ng^2*Ag/4)*(ig4/deltay4)^2;
G5=(mu0*Ng^2*Ag/4)*(ig5/deltay5)^2;
G6=(mu0*Ng^2*Ag/4)*(ig6/deltay6)^2;


%differential equations
ydot=vy;
vydot=(-G1+G2-G3+G4-G5+G6+m*g*sin(phi)-m*Ra*(wtheta^2)*cos(phi))/m;

%z - I'm using only this angle for two reason: because the other two are
%considered in the suspentions and also because this angles could not
%affect the distance from the guidance elector magnet if they're at the
%same height of the center of mass
psidot=wpsi;
wpsidot=(xB*(-G1+G2+G5-G6)+yB*(-P1+P2-P3+P4-P5+P6))/Ipsi;

% %equazioni elastiche
% qldot=vql;
% vqldot=-ql*w^2+gf*PHIf+gc*PHIc+gr*PHIr;

%LKT DEGLI ELETTROAMGNETI DELLA GUIDANCE
igdot1=ig1*deltaydot1/deltay1-(2/(mu0*Ng*alpha))*deltay1*(Ra*ig1-Vg1);
igdot2=ig2*deltaydot2/deltay2-(2/(mu0*Ng*alpha))*deltay2*(Ra*ig2-Vg2);
igdot3=ig3*deltaydot3/deltay3-(2/(mu0*Ng*alpha))*deltay3*(Ra*ig3-Vg3);
igdot4=ig4*deltaydot4/deltay4-(2/(mu0*Ng*alpha))*deltay4*(Ra*ig4-Vg4);
igdot5=ig5*deltaydot5/deltay5-(2/(mu0*Ng*alpha))*deltay5*(Ra*ig5-Vg5);
igdot6=ig6*deltaydot6/deltay6-(2/(mu0*Ng*alpha))*deltay6*(Ra*ig6-Vg6);


%diffetential equation
% ode = [ydot; vydot; psidot; wpsidot; qldot; vqldot; igdot1; igdot2; igdot3; igdot4; igdot5; igdot6];
ode = [ydot; vydot; psidot; wpsidot; igdot1; igdot2; igdot3; igdot4; igdot5; igdot6];

%compute the jacombian using the symbolic variables
a = jacobian(ode,stati);
b = jacobian(ode,input);

%set states linearization values
y      = 1;
vy     = rand(1);
psi    = rand(1);
wpsi   = rand(1);
% ql     = rand(1);
% vql    = rand(1);
ig1    = rand(1);
ig2    = rand(1);
ig3    = rand(1);
ig4    = rand(1);
ig5    = rand(1);
ig6    = rand(1);

%set inputs linearization values
Vg1 = rand(1);
Vg2 = rand(1);
Vg3 = rand(1);
Vg4 = rand(1);
Vg5 = rand(1);
Vg6 = rand(1);

%compute A B matricies
A=eval(a);
B=eval(b);

Co = ctrb(A,B); % controllability matrix
r = rank(Co);

% set Q and R matricies values respectively for states and input weights
weight_y      = 1;
weight_vy     = 0.0001;
weight_psi    = 0.0001;
weight_wpsi   = 0.0001;
% weight_ql     = 0.0001;
% weight_vql    = 0.0001;
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


K = lqr(A,B,Q,R);

%INITIAL VALUES

%set states inital values
y_init      = rand(1);
vy_init     = rand(1);
psi_init    = 0.01*rand(1);
wpsi_init   = 0.01*rand(1);
% ql_init     = rand(1);
% vql_init    = rand(1);
ig1_init    = 0;
ig2_init    = 0;
ig3_init    = 0;
ig4_init    = 0;
ig5_init    = 0;
ig6_init    = 0;



%-------REFERENCE SIGNAL GENERATION---------- DA FINIRE!!!!!!!!!!
ye      = 1;
vye     = 0;
psie    = 0;
wpsie   = 0;
qle     = 0;
vqle    = 0;
ig1e    = rand(1);
ig2e    = rand(1);
ig3e    = rand(1);
ig4e    = rand(1);
ig5e    = rand(1);
ig6e    = rand(1);

% xe = [ ye; vye; psie; wpsie; qle; vqle; ig1e; ig2e; ig3e; ig4e; ig5e; ig6e];
xe = [ ye; vye; psie; wpsie; ig1e; ig2e; ig3e; ig4e; ig5e; ig6e];

