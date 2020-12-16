%le parti mancanti sono segnate con i punti di sospensione
%le costanti non le ho ancora inserite
%gli array seguono l'ordine:
%1 Front Left 
%2 front right
%3 central l
%4 central r
%5 rear l
%6 rear r

%COSTANTI
zmax=1;
ymax=1;
PHI=1;
PHIf=1;
PHIc=1;
PHIr=1;
xB=1;
yB=1;
zB=1;
mu0=1;
Nl=1;
Ng=1;
Al=1;
Ag=1;
rho=1;
S=1;
Cd=1;
Cl=1;
m=1;
Iphi=1;
Itetha=1;
Ipsi=1;
g=9.81;
R=1;
l=1;
y0=1;
x0=1;

%STATO
x=state(1);
vx=state(2);
y=state(3);
vy=state(4);
z=state(5);
vz=state(6);
phi=state(7);
wphi=state(8);
tetha=state(9);
wtetha=state(10);
psi=state(11);
wpsi=state(12);
qv=state(13);
vqv=state(14);
ql=state(15);
vql=state(16);
il(1:6)=state(17:22);
ig(1:6)=state(23:28);
iL=state(29);
iR=state(30);

%VARIABILI NON APPARTENENTI ALLO STATO
deltaz=zeros(6,1);
deltay=zeros(6,1);
P=zeros(6,1);
L=zeros(6,1);
G=zeros(6,1);
state=zeros(20,1);

%spostamenti verticali dei rotori e dei sensori
deltaz(1)=zmax-z-phi*yB+tetha*xB-qv*PHI;
deltaz(2)=zmax-z+phi*yB+tetha*xB-qv*PHI;
deltaz(3)=zmax-z-phi*yB-qv*PHI;
deltaz(4)=zmax-z+phi*yB-qv*PHI;
deltaz(5)=zmax-z-phi*yB-tetha*xB-qv*PHI;
deltaz(6)=zmax-z+phi*yB-tetha*xB-qv*PHI;

%velocità degli spostamenti lungo l'asse z
deltazdot(1)=-vz-yB*wphi+xB*wtheta-vqv*PHI;
deltazdot(2)=-vz-yB*wphi-vqv*PHI;
deltazdot(3)=-vz-yB*wphi-xB*wtheta-vqv*PHI;
deltazdot(4)=-vz+yB*wphi+xB*wtheta-vqv*PHI;
deltazdot(5)=-vz+yB*wphi-vqv*PHI;
deltazdot(6)=-vz+yB*wphi-xB*wtheta-vqv*PHI;

%spostamenti orizzontali degli elettromagneti e dei sensori
deltay(2)=ymax-y-psi*xB-ql*PHI;
deltay(4)=ymax-y-ql*PHI;
deltay(6)=ymax-y+psi*xB-ql*PHI;
deltay(1)=2*ymax-deltay(2);
deltay(3)=2*ymax-deltay(4);
deltay(5)=2*ymax-deltay(6);

deltaydot(2)=-vy-xB*wpsi-vql*PHI;
deltaydot(4)=-vy-vql*PHI;
deltaydot(6)=-vy+xB*wpsi-vql*PHI;
deltaydot(1)=-deltaydot(2);
deltaydot(3)=-deltaydot(4);
deltaydot(5)=-deltaydot(6);

%forze di levitazione
L(1)=(mu0*Nl^2*Al/4)*(il(1)/deltaz(1))^2;
L(2)=(mu0*Nl^2*Al/4)*(il(2)/deltaz(2))^2;
L(3)=(mu0*Nl^2*Al/4)*(il(3)/deltaz(3))^2;
L(4)=(mu0*Nl^2*Al/4)*(il(4)/deltaz(4))^2;
L(5)=(mu0*Nl^2*Al/4)*(il(5)/deltaz(5))^2;
L(6)=(mu0*Nl^2*Al/4)*(il(6)/deltaz(6))^2;

%forse della guidance
G(1)=(mu0*Ng^2*Ag/4)*(ig(1)/deltay(1))^2;
G(2)=(mu0*Ng^2*Ag/4)*(ig(2)/deltay(2))^2;
G(3)=(mu0*Ng^2*Ag/4)*(ig(3)/deltay(3))^2;
G(4)=(mu0*Ng^2*Ag/4)*(ig(4)/deltay(4))^2;
G(5)=(mu0*Ng^2*Ag/4)*(ig(5)/deltay(5))^2;
G(6)=(mu0*Ng^2*Ag/4)*(ig(6)/deltay(6))^2;

%propulsione
P(1:6)=[0 0 0 0 0 0];.......

%distribuzioni delle forze sul modello elastico
lf=(L(1)+L(2)-((rho*S*Cl/2)*vx^2+m*R*(wtetha^2)*sin(phi)+m*g*cos(phi))/3)/l;
lc=(L(3)+L(4)-((rho*S*Cl/2)*vx^2+m*R*(wtetha^2)*sin(phi)+m*g*cos(phi))/3)/l;
lr=(L(5)+L(6)-((rho*S*Cl/2)*vx^2+m*R*(wtetha^2)*sin(phi)+m*g*cos(phi))/3)/l;
gf=(G(1)-G(2)+(m*g*sin(phi)-m*R*cos(phi)*wtetha^2)/3)/l;
gc=(G(3)-G(4)+(m*g*sin(phi)-m*R*cos(phi)*wtetha^2)/3)/l;
gr=(G(5)-G(6)+(m*g*sin(phi)-m*R*cos(phi)*wtetha^2)/3)/l;

%EQUAZIONI CARDINALI DELLA MECCANICA
xdot=vx;
xdotdot=(P(1)+P(2)+P(3)+P(4)+P(5)+P(6)-(rho*S*Cd/2)*vx^2)/m;
ydot=vy;
ydotdot=(-G(1)+G(2)-G(3)+G(4)-G(5)+G(6)+m*g*sin(phi)-m*R*(wtetha^2)*cos(phi))/m;
zdot=vz;
zdotdot=(L(1)+L(2)+L(3)+L(4)+L(5)+L(6)-(rho*S*Cl/2)*vx^2-m*R*(wtetha^2)*sin(phi)-m*g*cos(phi))/m;
phidot=wphi;
phidotdot=(yB*(L(1)-L(2)+L(3)-L(4)+L(5)-L(6))+zB*(G(1)-G(2)+G(3)-G(4)+G(5)-G(6)))/Iphi;
tethadot=wtetha;
tethadotdot=(xB*(-L(1)-L(2)+L(5)+L(6))+zB*(P(1)+P(2)+P(3)+P(4)+P(5)+P(6))-zA*(rho*S*Cd/2)*vx^2+xA*(rho*S*Cl/2)*vx^2)/Itetha;
psidot=wpsi;
psidotdot=(xB*(-G(1)+G(2)+G(5)-G(6))+yB*(-P(1)+P(2)-P(3)+P(4)-P(5)+P(6)))/Ipsi;

%EQUAZIONI DEL MODELLO ELASTICO (per ora ne ho scritte solo due, poi vediamo)
qvdot=vqv;
qvdotdot=-qv*w^2+lf*PHIf+lc*PHIc+lr*PHIr;
qldot=vql;
qldotdot=-qlw^2+gf*PHIf+gc*PHIc+gr*PHIr;

%LKT DEI ROTORI
ildot(1:6)=il(1:6).*;.......

ildot(1)=il(1)*
    
%LKT DEGLI ELETTROAMGNETI DELLA GUIDANCE
igdot(1:6)=;........
igdot(1)=
%LKT DEGLI STATORI
iRdot=0;........
iLdot=0;........