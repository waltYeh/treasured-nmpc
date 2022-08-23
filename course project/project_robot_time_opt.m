function project_robot_time_opt()
addpath('../../casadi')
import casadi.*;

clear;
clc;

x0=[-5;-4;0;0];
xF=[pi/2;0;0;0];
% Q=diag([1,1,1,1]);
R=diag([0,1e-6]);
% Q=diag([100,100,10,10]);
% R=diag([1000,10]);

[~, ~, ~, ~, state_constr, input_constr] = project_parameters;
N=100;
opti=casadi.Opti();
X=opti.variable(4,N+1);
U=opti.variable(2,N);
tf=opti.variable(1,1);
tau=0;
J=0;
for i=1:N
    tau=tau+1/N;
%     x_next = rk4(@(t,x,u)diff_eq(t,x,u),tf/N,tau,X(:,i),U(:,i));
    x_next = rk4(@(t,x,u)diff_eq_time_trans(t,x,u,tf),1/N,tau,X(:,i),U(:,i));
    opti.subject_to(X(:,i+1)==x_next);
    J = J + U(:,i)' * R * U(:,i)*1/N;
end
opti.subject_to(X(:,1) == x0);
opti.subject_to(X(:,end) == xF);
opti.subject_to(-state_constr <= X(3,:) <= state_constr);
opti.subject_to(-state_constr <= X(4,:) <= state_constr);
opti.subject_to(-input_constr <= U(1,:) <= input_constr);
opti.subject_to(-input_constr <= U(2,:) <= input_constr);
opti.subject_to(tf > 0);
opti.minimize(tf+J);
opti.solver('ipopt');
opti.set_initial(tf,3);
sol = opti.solve();
q1=sol.value(X(1,:))';
q2=sol.value(X(2,:))';
w1=sol.value(X(3,:))';
w2=sol.value(X(4,:))';
u1=sol.value(U(1,:))';
u2=sol.value(U(2,:))';
tx=linspace(0,sol.value(tf(:)),N+1);
tu=linspace(0,sol.value(tf(:))*N/(N+1),N);
figure(1)
subplot(3,1,1)
plot(tx,q1,tx,q2)
title('states')
xlabel('t (s)')
ylabel('angle (rad)')
legend('q_1','q_2')
subplot(3,1,2)
plot(tx,w1,tx,w2)
title('states')
xlabel('t (s)')
ylabel('angular velocity (rad/s)')
legend('omiga_1','omiga_2')
subplot(3,1,3)
plot(tu,u1,tu,u2)
title('inputs')
xlabel('t (s)')
ylabel('torque (Nm)')
legend('u_1','u_2')

end
function xf = rk4(ode,h,t,x,u)
  k1 = ode(t,x,u);
  k2 = ode(t,x+h/2*k1,u);
  k3 = ode(t,x+h/2*k2,u);
  k4 = ode(t,x+h*k3,  u);
  xf = x + h/6 * (k1 + 2*k2 + 2*k3 + k4); 
end
function dx=diff_eq_time_trans(tau,x,u,tf)
q1=x(1);
q2=x(2);
w1=x(3);
w2=x(4);
[b,c,g,~,~,~]=project_parameters;
B=[b(1)+b(2)*cos(q2), b(3)+b(4)*cos(q2);b(3)+b(4)*cos(q2),b(5)];
C=-c*sin(q2)*[w1,w1+w2;-w1,0];
G=[g(1)*cos(q1)+g(2)*cos(q1+q2);g(2)*cos(q1+q2)];
dq1=w1;
dq2=w2;
dw=B\(u-G-C*[w1;w2]);
dw1=dw(1);
dw2=dw(2);
dx=tf*[dq1;dq2;dw1;dw2];
end
function dx=diff_eq(t,x,u)
q1=x(1);
q2=x(2);
w1=x(3);
w2=x(4);
[b,c,g,~,~,~]=project_parameters;
B=[b(1)+b(2)*cos(q2), b(3)+b(4)*cos(q2);b(3)+b(4)*cos(q2),b(5)];
C=-c*sin(q2)*[w1,w1+w2;-w1,0];
G=[g(1)*cos(q1)+g(2)*cos(q1+q2);g(2)*cos(q1+q2)];
dq1=w1;
dq2=w2;
dw=B\(u-G-C*[w1;w2]);
dw1=dw(1);
dw2=dw(2);
dx=[dq1;dq2;dw1;dw2];
end
function [b, c, g, l, state_constr, input_constr] = project_parameters
%% Definition of system parameters
b=[200;50;23.5;25;122.5];
c=-25;
g=[784.8;245.3];
l=[0.5;0.5];

%% Constraints
state_constr=3/2*pi;
input_constr=1000;

%% Scaling factors

end

