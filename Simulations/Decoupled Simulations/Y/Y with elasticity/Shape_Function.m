function [shape,int11,int21,int31,omega_1]=Shape_Function()
clc
close all
clear
%%% Decommentare i grafici se si vogliono visualizzare
syms  k
l=8;    %lunghezza singolo vagone in metri 
%imposizione condizioni al contorno
a=[1, 1, -1, 0;
    1, -1, 0, -1; 
    exp(k*3*l), exp(-k*3*l), -cos(k*3*l), -sin(k*3*l);
    exp(k*3*l), -exp(-k*3*l), sin(k*3*l), -cos(k*3*l)];
D=simplify(det(a));

% calcolo K1 da det(A)=0
% figure('Name','Zeri del determinante')  
% fplot([D 0]);
% xlim([0,1]);
% ylim([-100,100]);
K1=vpasolve(D==0,k,0.2);


%calcolo frequenze naturali dei 3 modi
rho= 7833; %densità acciaio in Kg/m^3
E=200*10^9; %modulo di Young acciaio in GPa
A=0.3; %sezione trave in m^2
I=0.44; %inerzia orizzontale in m^4
omega_1=sqrt((E*I)/(rho*A))*K1^2;


psi=(sin(3*K1*l)+cos(3*K1*l)-exp(3*K1*l))/(sin(3*K1*l)-cos(3*K1*l)+exp(3*K1*l)); 


%esprimo tutto in funzione di F e impongo la condizione di ortogonalità
%per determinarlo

%%%%%%%%% MODO 1 %%%%%%
syms s
expr1= ((psi+1)*exp(K1*s)+(psi-1)*exp(-K1*s)+2*psi*cos(K1*s)+2*sin(K1*s))^2;
integrale1=vpa(int(expr1,s,0,3*l));
F1=sqrt(4/(rho*A*integrale1));   
E1=F1*psi;   
D1=(E1-F1)/2;  
C1=(E1+F1)/2;  

shape=C1*exp(K1*s)+D1*exp(-K1*s)+E1*cos(K1*s)+F1*sin(K1*s);
% figure('Name','Shape function 1')
% fplot([shape 0]);
% xlim([0,24]);
% ylim([-10^-1,10^-1]);

 
% %%%%%% SOLUZIONE PARTICOLARE %%%%%%
% %MODO 1
int11=int(shape/l,s,0,l);
int21=int(shape/l,s,l,2*l);
int31=int(shape/l,s,2*l,3*l);
