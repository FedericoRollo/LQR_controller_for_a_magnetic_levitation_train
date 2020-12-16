
%------------------------STATES STATEMENT--------------------------------

syms x vx y vy z vz phi wphi theta wtheta psi wpsi qv vqv ql vql ir1 ir2 ir3 ir4 ir5 ir6 ig1 ig2 ig3 ig4 ig5 ig6 iqL idL iqR idR
stati = [ x vx y vy z vz phi wphi theta wtheta psi wpsi qv vqv ql vql ir1 ir2 ir3 ir4 ir5 ir6 ig1 ig2 ig3 ig4 ig5 ig6 iqL idL iqR idR];

%-----------------------INPUT STATEMENT----------------------------------

syms Vr1 Vr2 Vr3 Vr4 Vr5 Vr6 Vg1 Vg2 Vg3 Vg4 Vg5 Vg6 Vql Vdl Vqr Vdr
input = [Vr1 Vr2 Vr3 Vr4 Vr5 Vr6 Vg1 Vg2 Vg3 Vg4 Vg5 Vg6 Vql Vdl Vqr Vdr];

%gli array seguono l'ordine:
%1 front left 
%2 front right
%3 central left
%4 central right
%5 rear left
%6 rear right

%-----------------------------CONSTANTS-------------------------------------
w=1;
PHI=1;
PHIf=1;
PHIc=1;
PHIr=1;
xB=1;
yB=1;
zB=1;
xA=1;
zA=1;
mu0=1;
Nr=1;
Ng=1;
Ns=1;
tau=1;
alpha=1; 
Al=1;
Ag=1;
rho=1;
S=1;
Cd=1;
Cl=1;
m=1;
Iphi=1;
Itheta=1;
Ipsi=1;
g=9.81;
R=1;
l=1;

%--------------------SUPPORT ALGEBRIC FUNCTIONS------------------------

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

%spostamenti orizzontali degli elettromagneti e dei sensori
deltay2=-y-psi*xB-ql*PHI;
deltay4=-y-ql*PHI;
deltay6=-y+psi*xB-ql*PHI;
deltay1=-deltay2;
deltay3=-deltay4;
deltay5=-deltay6;

deltaydot2=-vy-wpsi*xB-vql*PHI;
deltaydot4=-vy-vql*PHI;
deltaydot6=-vy+wpsi*xB-vql*PHI;
deltaydot1=-deltaydot2;
deltaydot3=-deltaydot4;
deltaydot5=-deltaydot6;

%forze di levitazione
L1=(mu0*Nr^2*Al/4)*(ir1/deltaz1)^2;
L2=(mu0*Nr^2*Al/4)*(ir2/deltaz2)^2;
L3=(mu0*Nr^2*Al/4)*(ir3/deltaz3)^2;
L4=(mu0*Nr^2*Al/4)*(ir4/deltaz4)^2;
L5=(mu0*Nr^2*Al/4)*(ir5/deltaz5)^2;
L6=(mu0*Nr^2*Al/4)*(ir6/deltaz6)^2;

%forse della guidance
G1=(mu0*Ng^2*Ag/4)*(ig1/deltay1)^2;
G2=(mu0*Ng^2*Ag/4)*(ig2/deltay2)^2;
G3=(mu0*Ng^2*Ag/4)*(ig3/deltay3)^2;
G4=(mu0*Ng^2*Ag/4)*(ig4/deltay4)^2;
G5=(mu0*Ng^2*Ag/4)*(ig5/deltay5)^2;
G6=(mu0*Ng^2*Ag/4)*(ig6/deltay6)^2;

%propulsione
P1=(((27*pi)*(mu0*Ns*Nr*alpha)^2)/(32*tau))*(ir1*iqL/(deltaz1^2));
P2=(((27*pi)*(mu0*Ns*Nr*alpha)^2)/(32*tau))*(ir2*iqR/(deltaz2^2));
P3=(((27*pi)*(mu0*Ns*Nr*alpha)^2)/(32*tau))*(ir3*iqL/(deltaz3^2));
P4=(((27*pi)*(mu0*Ns*Nr*alpha)^2)/(32*tau))*(ir4*iqR/(deltaz4^2));
P5=(((27*pi)*(mu0*Ns*Nr*alpha)^2)/(32*tau))*(ir5*iqL/(deltaz5^2));
P6=(((27*pi)*(mu0*Ns*Nr*alpha)^2)/(32*tau))*(ir6*iqR/(deltaz6^2));

%distribuzioni delle forze sul modello elastico
lf=(L1+L2-((rho*S*Cl/2)*(vx)^2+m*R*(wtheta^2)*sin(phi)+m*g*cos(phi))/3)/l;
lc=(L3+L4-((rho*S*Cl/2)*(vx)^2+m*R*(wtheta^2)*sin(phi)+m*g*cos(phi))/3)/l;
lr=(L5+L6-((rho*S*Cl/2)*(vx)^2+m*R*(wtheta^2)*sin(phi)+m*g*cos(phi))/3)/l;
gf=(G1-G2+(m*g*sin(phi)-m*R*cos(phi)*wtheta^2)/3)/l;
gc=(G3-G4+(m*g*sin(phi)-m*R*cos(phi)*wtheta^2)/3)/l;
gr=(G5-G6+(m*g*sin(phi)-m*R*cos(phi)*wtheta^2)/3)/l;

%----------------------DIFFERENTIAL EQUATIONS------------------------------

%EQUAZIONI CARDINALI DELLA MECCANICA
xdot=(vx);
vxdot=(P1+P2+P3+P4+P5+P6-(rho*S*Cd/2)*(vx)^2)/m;

ydot=vy;
vydot=(-G1+G2-G3+G4-G5+G6+m*g*sin(phi)-m*R*(wtheta^2)*cos(phi))/m;

zdot=vz;
vzdot=(L1+L2+L3+L4+L5+L6-(rho*S*Cl/2)*(vx)^2-m*R*(wtheta^2)*sin(phi)-m*g*cos(phi))/m;
%x
phidot=wphi;
wphidot=(yB*(L1-L2+L3-L4+L5-L6)+zB*(G1-G2+G3-G4+G5-G6))/Iphi;
%y
thetadot=wtheta;
wthetadot=(xB*(-L1-L2+L5+L6)+zB*(P1+P2+P3+P4+P5+P6)-zA*(rho*S*Cd/2)*(vx)^2+xA*(rho*S*Cl/2)*(vx)^2)/Itheta;
%z
psidot=wpsi;
wpsidot=(xB*(-G1+G2+G5-G6)+yB*(-P1+P2-P3+P4-P5+P6))/Ipsi;

%EQUAZIONI DEL MODELLO ELASTICO (per ora ne ho scritte solo due, poi vediamo)
qvdot=vqv;
vqvdot=-qv*w^2+lf*PHIf+lc*PHIc+lr*PHIr;
qldot=vql;
vqldot=-ql*w^2+gf*PHIf+gc*PHIc+gr*PHIr;

%LKT DEI ROTORI
irdot1=ir1*deltazdot1/deltaz1-(2/(mu0*Nr*alpha))*deltaz1*(R*ir1-Vr1);
irdot2=ir2*deltazdot2/deltaz2-(2/(mu0*Nr*alpha))*deltaz2*(R*ir2-Vr2);
irdot3=ir3*deltazdot3/deltaz3-(2/(mu0*Nr*alpha))*deltaz3*(R*ir3-Vr3);
irdot4=ir4*deltazdot4/deltaz4-(2/(mu0*Nr*alpha))*deltaz4*(R*ir4-Vr4);
irdot5=ir5*deltazdot5/deltaz5-(2/(mu0*Nr*alpha))*deltaz5*(R*ir5-Vr5);
irdot6=ir6*deltazdot6/deltaz6-(2/(mu0*Nr*alpha))*deltaz6*(R*ir6-Vr6);
    
%LKT DEGLI ELETTROAMGNETI DELLA GUIDANCE
igdot1=ig1*deltaydot1/deltay1-(2/(mu0*Ng*alpha))*deltay1*(R*ig1-Vg1);
igdot2=ig2*deltaydot2/deltay2-(2/(mu0*Ng*alpha))*deltay2*(R*ig2-Vg2);
igdot3=ig3*deltaydot3/deltay3-(2/(mu0*Ng*alpha))*deltay3*(R*ig3-Vg3);
igdot4=ig4*deltaydot4/deltay4-(2/(mu0*Ng*alpha))*deltay4*(R*ig4-Vg4);
igdot5=ig5*deltaydot5/deltay5-(2/(mu0*Ng*alpha))*deltay5*(R*ig5-Vg5);
igdot6=ig6*deltaydot6/deltay6-(2/(mu0*Ng*alpha))*deltay6*(R*ig6-Vg6);

%LKT DEGLI STATORI
iqLdot=iqL*(deltazdot1/deltaz1+deltazdot3/deltaz3+deltazdot5/deltaz5)-((3*pi*mu0*(Nr^2)*alpha)/(4*tau))*(vx)*(ir1/deltaz1+ir3/deltaz3+ir5/deltaz5)+(4/(3*mu0*(Ns^2)*alpha))*(deltaz1+deltaz3+deltaz5)*(Vql-R*iqL);
iqRdot=iqR*(deltazdot2/deltaz2+deltazdot4/deltaz4+deltazdot6/deltaz6)-((3*pi*mu0*(Nr^2)*alpha)/(4*tau))*(vx)*(ir2/deltaz2+ir4/deltaz4+ir6/deltaz6)+(4/(3*mu0*(Ns^2)*alpha))*(deltaz2+deltaz4+deltaz6)*(Vqr-R*iqR);
idLdot=idL*(deltazdot1/deltaz1+deltazdot3/deltaz3+deltazdot5/deltaz5)-(pi/tau)*(vx)*(iqL)+((3*pi*mu0*(Nr^2)*alpha)/(4*tau))*(2*(ir1*deltazdot1/deltaz1+ir3*deltazdot3/deltaz3+ir5*deltazdot5/deltaz5)-(irdot1/deltaz1+irdot3/deltaz3+irdot5/deltaz5))+(4/(3*mu0*(Ns^2)*alpha))*(deltaz1+deltaz3+deltaz5)*(Vdl-R*idL);
idRdot=idR*(deltazdot2/deltaz2+deltazdot4/deltaz4+deltazdot6/deltaz6)-(pi/tau)*(vx)*(iqR)+((3*pi*mu0*(Nr^2)*alpha)/(4*tau))*(2*(ir2*deltazdot2/deltaz2+ir4*deltazdot4/deltaz4+ir6*deltazdot6/deltaz6)-(irdot2/deltaz2+irdot4/deltaz4+irdot6/deltaz6))+(4/(3*mu0*(Ns^2)*alpha))*(deltaz2+deltaz4+deltaz6)*(Vdr-R*idR);

%DIFFERENTIAL EQUATION VECTOR
ode = [xdot; vxdot; ydot; vydot; zdot; vzdot; phidot; wphidot; thetadot; wthetadot; psidot; wpsidot; qvdot; vqvdot; qldot; vqldot; irdot1; irdot2; irdot3; irdot4; irdot5; irdot6; igdot1; igdot2; igdot3; igdot4; igdot5; igdot6; iqLdot; idLdot; iqRdot; idRdot];

%----------------------------LINEARIZATION---------------------------------

%compute the jacombian using the symbolic variables
a = jacobian(ode,stati);
b = jacobian(ode,input);

%set states linearization values
x      = rand(1);
vx     = 200;
y      = 1;
vy     = rand(1);
z      = 1;
vz     = rand(1);
phi    = rand(1);
wphi   = rand(1);
theta  = rand(1);
wtheta = rand(1);
psi    = rand(1);
wpsi   = rand(1);
qv     = rand(1);
vqv    = rand(1);
ql     = rand(1);
vql    = rand(1);
ir1    = rand(1);
ir2    = rand(1);
ir3    = rand(1);
ir4    = rand(1);
ir5    = rand(1);
ir6    = rand(1);
ig1    = rand(1);
ig2    = rand(1);
ig3    = rand(1);
ig4    = rand(1);
ig5    = rand(1);
ig6    = rand(1);
iqL    = rand(1);
idL    = rand(1);
iqR    = rand(1);
idR    = rand(1);

%set inputs linearization values
Vr1 = rand(1);
Vr2 = rand(1);
Vr3 = rand(1);
Vr4 = rand(1);
Vr5 = rand(1);
Vr6 = rand(1);
Vg1 = rand(1);
Vg2 = rand(1);
Vg3 = rand(1);
Vg4 = rand(1);
Vg5 = rand(1);
Vg6 = rand(1);
Vql = rand(1);
Vdl = rand(1);
Vqr = rand(1);
Vdr = rand(1);

%compute A B matricies
A=eval(a);
B=eval(b);

Co = ctrb(A,B); % controllability matrix
r = rank(Co);

%--------------------------LQR CONTROL------------------------------------

% set Q and R matricies values respectively for states and input weights
weight_x      = 0.0001;
weight_vx     = 100;
weight_y      = 1;
weight_vy     = 0.0001;
weight_z      = 100;
weight_vz     = 0.0001;
weight_phi    = 0.0001;
weight_wphi   = 0.0001;
weight_theta  = 0.0001;
weight_wtheta = 0.0001;
weight_psi    = 0.0001;
weight_wpsi   = 0.0001;
weight_qv     = 0.0001;
weight_vqv    = 0.0001;
weight_ql     = 0.0001;
weight_vql    = 0.0001;
weight_ir1    = 0.0001;
weight_ir2    = 0.0001;
weight_ir3    = 0.0001;
weight_ir4    = 0.0001;
weight_ir5    = 0.0001;
weight_ir6    = 0.0001;
weight_ig1    = 0.0001;
weight_ig2    = 0.0001;
weight_ig3    = 0.0001;
weight_ig4    = 0.0001;
weight_ig5    = 0.0001;
weight_ig6    = 0.0001;
weight_iqL    = 0.0001;
weight_idL    = 0.0001;
weight_iqR    = 0.0001;
weight_idR    = 0.0001;

Q = diag([weight_x weight_vx weight_y weight_vy weight_z weight_vz weight_phi weight_wphi weight_theta weight_wtheta weight_psi weight_wpsi weight_qv weight_vqv weight_ql weight_vql weight_ir1 weight_ir2 weight_ir3 weight_ir4 weight_ir5 weight_ir6 weight_ig1 weight_ig2 weight_ig3 weight_ig4 weight_ig5 weight_ig6 weight_iqL weight_idL weight_iqR weight_idR]);

weight_Vr1 = 0.0001;
weight_Vr2 = 0.0001;
weight_Vr3 = 0.0001;
weight_Vr4 = 0.0001;
weight_Vr5 = 0.0001;
weight_Vr6 = 0.0001;
weight_Vg1 = 0.0001;
weight_Vg2 = 0.0001;
weight_Vg3 = 0.0001;
weight_Vg4 = 0.0001;
weight_Vg5 = 0.0001;
weight_Vg6 = 0.0001;
weight_Vql = 0.0001;
weight_Vdl = 0.0001;
weight_Vqr = 0.0001;
weight_Vdr = 0.0001;

R = diag([weight_Vr1 weight_Vr2 weight_Vr3 weight_Vr4 weight_Vr5 weight_Vr6 weight_Vg1 weight_Vg2 weight_Vg3 weight_Vg4 weight_Vg5 weight_Vg6 weight_Vql weight_Vdl weight_Vqr weight_Vdr]);


K = lqr(A,B,Q,R);

%-------REFERENCE SIGNAL GENERATION----------
xe      = 0;
vxe     = 100/3.6;
ye      = 1;
vye     = 0;
ze      = 1;
vze     = 0;
phie    = 0;
wphie   = 0;
thetae  = 0;
wthetae = 0;
psie    = 0;
wpsie   = 0;
qve     = 0;
vqve    = 0;
qle     = 0;
vqle    = 0;
ir1e    = 0;
ir2e    = 0;
ir3e    = 0;
ir4e    = 0;
ir5e    = 0;
ir6e    = 0;
ig1e    = 0;
ig2e    = 0;
ig3e    = 0;
ig4e    = 0;
ig5e    = 0;
ig6e    = 0;
iqLe    = 0;
idLe    = 0;
iqRe    = 0;
idRe    = 0;

xe = [xe; vxe; ye; vye; ze; vze; phie; wphie; thetae; wthetae; psie; wpsie; qve; vqve; qle; vqle; ir1e; ir2e; ir3e; ir4e; ir5e; ir6e; ig1e; ig2e; ig3e; ig4e; ig5e; ig6e; iqLe; idLe; iqRe; idRe];
    

%output=[xdot, xdotdot, ydot, ydotdot, zdot, zdotdot, phidot, phidotdot, thetadot, thetadotdot, psidot, psidotdot, qvdot, qvdotdot, qldot, qldotdot, irdot', igdot', iqLdot, iqRdot, idLdot, idRdot]';