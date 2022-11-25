clear all
clc

%% system parameters
M = 0.881; m = 0.015; R = 0.025; D = 0.137; L = 0.040;
Jw = 4.6875*10^(-6); Jp = 4.698667*10^(-4); g=9.8;

Kt = 0.02; Kb=0.02; Ra=1;
Kr = Kt/Ra;

syms Vr Vl theta theta_d x x_d delta delta_d 

%% solve Ax=B
A = [0.5*M*R+(m*R^2+Jw)/R , 0.5*R*L*M*cos(theta) , (D*(m*R^2+Jw))/(2*R);
     0.5*M*R+(m*R^2+Jw)/R , 0.5*R*L*M*cos(theta) , (-D*(m*R^2+Jw))/(2*R);
     2*L*M*cos(theta) , M*L^2*(cos(theta))^2 + Jp , 0];

B = [Kr*Vl-Kr*Kb*theta_d-0.5*R*L*M*theta_d^2*sin(theta);
     Kr*Vr-Kr*Kb*theta_d-0.5*R*L*M*theta_d^2*sin(theta);
     L*M*theta_d*sin(theta)*x_d+M*g*L*sin(theta)+M*L^2*theta_d^2*sin(theta)*cos(theta)-Kr*Vr-Kr*Vl+Kr*Kb*theta_d];

X = inv(A) * B;
X = simplify(X);


%% linearization
states = [x  theta delta  x_d theta_d delta_d];
f = [x_d;theta_d;delta_d; X]';

%A
As=jacobian(f,states); 
A=subs(As , [[states], [Vr, Vl]] , [0 0 0 0 0 0 0 0]);
A = double(A);

%B
Bs=jacobian([f],[Vr, Vl]); 
B=subs(Bs , [[states], [Vr ,Vl]] , [0 0 0 0 0 0 0 0]);
B = double(B);

%C
C=[1 0 0 0 0 0;
   0 1 0 0 0 0;
   0 0 1 0 0 0;
   0 0 0 1 0 0;
   0 0 0 0 1 0;
   0 0 0 0 0 1];




%% controllability  &  visibility  &    instability

% controllability
M1=rank(ctrb(A,B))

% visibility
M2=rank(obsv(A,C))

%  instability
E=eig(A)