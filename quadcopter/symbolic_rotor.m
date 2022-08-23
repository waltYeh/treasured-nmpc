function [A_num,B_num,B_double_thrust,C,D]=symbolic_rotor()

syms mg mq L betta Ig Iq f theta M2 tau t gravity theta_dot beta_dot x_dot z_dot x z

D=[mg+mq,0,0,-mg*L*sin(betta);
    0,mg+mq,0,-mg*L*cos(betta);
    0,0,Iq,0;
    -mg*L*sin(betta),-mg*L*cos(betta),0,mg*L*L+Ig];
F=[f*sin(theta);f*cos(theta);M2-tau;tau];
C=[0,0,0,-L*mg*cos(betta)*beta_dot;0,0,0,L*mg*sin(betta)*beta_dot;0,0,0,0;0,0,0,0];
G=[0;(mg+mq)*gravity;0;-mg*gravity*L*cos(betta)];
q_dot = [x_dot; z_dot;theta_dot; beta_dot];
f5678=D\(F-C*q_dot-G);
f5=f5678(1);
f6=f5678(2);
f7=f5678(3);
f8=f5678(4);
JA = [0,0,0,0,1,0,0,0;
    0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,1;
    diff(f5,x),diff(f5,z),diff(f5,theta),diff(f5,betta),diff(f5,x_dot),diff(f5,z_dot),diff(f5,theta_dot),diff(f5,beta_dot);
    diff(f6,x),diff(f6,z),diff(f6,theta),diff(f6,betta),diff(f6,x_dot),diff(f6,z_dot),diff(f6,theta_dot),diff(f6,beta_dot);
    diff(f7,x),diff(f7,z),diff(f7,theta),diff(f7,betta),diff(f7,x_dot),diff(f7,z_dot),diff(f7,theta_dot),diff(f7,beta_dot);
    diff(f8,x),diff(f8,z),diff(f8,theta),diff(f8,betta),diff(f8,x_dot),diff(f8,z_dot),diff(f8,theta_dot),diff(f8,beta_dot)];
x0=zeros(8,1);
x0(4)=pi/2;
x_sym=[x;z;theta;betta;x_dot;z_dot;theta_dot;beta_dot];
u_sym = [f;M2;tau];
u0=[(gravity*(mg+mq));0;0];
A=subs(JA,x_sym,x0);
A=subs(A,u_sym,u0);
A=simplify(A);
A=simplify(A);
JB=[0,0,0;
    0,0,0;
    0,0,0;
    0,0,0;
    diff(f5,f),diff(f5,M2),diff(f5,tau);
    diff(f6,f),diff(f6,M2),diff(f6,tau);
    diff(f7,f),diff(f7,M2),diff(f7,tau);
    diff(f8,f),diff(f8,M2),diff(f8,tau)];
B=subs(JB,x_sym,x0);
B=subs(B,u_sym,u0);
B=simplify(B);
B=simplify(B);
C=[1,0,0,0, 0,0,0,0;
    0,1,0,0, 0,0,0,0;
    0,0,0,1, 0,0,0,0];
[mg_num, mq_num, Ig_num, Iq_num, L_center_of_mass, ~, gravity_num] = project_parameters(0);
const_num = [mg_num, mq_num, Ig_num, Iq_num, L_center_of_mass, gravity_num];
const_sym = [mg, mq, Ig, Iq, L, gravity];
A_num=double(subs(A,const_sym,const_num));
B_num=double(subs(B,const_sym,const_num));
B_double_thrust = B_num;
B_double_thrust(6,1) = 2*B_double_thrust(6,1);
B_double_thrust(6,2) = B_double_thrust(6,1);
B_double_thrust(7,1) = -B_double_thrust(7,2);
D=zeros(3,3);
G=tf(ss(A_num,B_num,C,D));
% zero(G(1,2))
% zero(G(1,3))
% zero(G(3,3))
% G=ss2tf(A,B,C,D)
% G=C*inv(s*eye(8)-A)*B;
G_reform = [G(1,2),0,G(1,3);0,G(2,1),0;G(3,2),0,G(3,3)];
G_D = [G_reform(1,1),0,0;0,G_reform(2,2),0;0,0,G_reform(3,3)];
syms GK1 GK2 GK3
GK=diag([GK1,GK2,GK3]);
% F=G_D*GK;
% R=G_reform\F;
assignin('base','A',A_num);
assignin('base','B',B_num);
assignin('base','B_dbth',B_double_thrust);
assignin('base','C',C);
assignin('base','D',D);
end
