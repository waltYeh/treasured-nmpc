function project_robot_opt()
addpath('../../casadi')
import casadi.*;
clear;
clc;
close
x0=[-5;-4;0;0];
% x0=[0;0;0;0];
xF=[pi/2;0;0;0];
% Q=diag([1,1,1,1]);
% R=diag([1,1]);
Q=diag([0,0,0,0]);
R=diag([10,10]);

[~, ~, ~, arm_len, state_constr, input_constr] = project_parameters;
N=30;
t0=0;
tF=3;
Ts=(tF-t0)/N;
opti=casadi.Opti();
X=opti.variable(4,N+1);
U=opti.variable(2,N);
W=5e8;
J=0;
for i=1:N
    x_next = rk4(@(t,x,u)diff_eq(t,x,u),Ts,0,X(:,i),U(:,i));
    opti.subject_to(X(:,i+1)==x_next);
    error = x_next - xF;
    J = J + (U(:,i)' * R * U(:,i)+error'*Q*error)*Ts;
end
opti.subject_to(X(:,1) == x0);
if W==0
    opti.subject_to(X(:,end) == xF);
end
J = J + (X(:,end) - xF)'*(X(:,end) - xF)*W;
% opti.set_initial(W,1e+7);
opti.subject_to(-state_constr <= X(3,:) <= state_constr);
opti.subject_to(-state_constr <= X(4,:) <= state_constr);
opti.subject_to(-input_constr <= U(1,:) <= input_constr);
opti.subject_to(-input_constr <= U(2,:) <= input_constr);

% free_zone=round(tF/12/Ts);
% i=round(tF/6/Ts);
% opti.subject_to(X(3,free_zone:i)>=X(4,free_zone:i));
% i=N-round(tF/6/Ts);
% opti.subject_to(X(3,i:end-free_zone)<=X(4,i:end-free_zone));


% ocp.set_initial(X,
opti.minimize(J);
opti.solver('ipopt');
sol = opti.solve();
q1=sol.value(X(1,:))';
q2=sol.value(X(2,:))';
w1=sol.value(X(3,:))';
w2=sol.value(X(4,:))';
u1=sol.value(U(1,:))';
u2=sol.value(U(2,:))';

unweighted_terminal_cost=([q1(end);q2(end);w1(end);w2(end)] - xF)'*([q1(end);q2(end);w1(end);w2(end)] - xF);
energy_cost=sol.value(J)-unweighted_terminal_cost*W;
tx=linspace(t0,tF,N+1);
tu=linspace(t0,tF-Ts,N)';
figure(2)
subplot(3,1,1)
plot(tx,q1,tx,q2)
axis([0 3 -5 5])
title('states')
xlabel('t (s)')
ylabel('angle (rad)')
legend('q_1','q_2')

subplot(3,1,2)
plot(tx,w1,tx,w2)
axis([0 3 -5 5])
title('states')
xlabel('t (s)')
ylabel('angular velocity (rad/s)')
legend('omiga_1','omiga_2')

subplot(3,1,3)
stairs([tu,tu],[u1,u2])
axis([0 3 -1000 1000])
% hold on
% stairs(tu,u2,'r')
title('inputs')
xlabel('t (s)')
ylabel('torque (Nm)')
legend('u_1','u_2')
figure(3)
subplot(2,1,1)
plot(q1,q2,'-+')
x1=arm_len(1)*cos(q1);
y1=arm_len(1)*sin(q1);
x2=x1+arm_len(2)*cos(q1+q2);
y2=y1+arm_len(2)*sin(q1+q2);
axis equal
title('q_1-q_2-plane')
xlabel('q_1(rad)')
ylabel('q_2(rad)')
axis([-6,2,-5,0])
subplot(2,1,2)
% plot(x1,y1,'b--*',x2,y2,'r--*')
hold on
for i=1:length(q1)
    plot([0,x1(i)],[0,y1(i)],'b')
    plot([x1(i),x2(i)],[y1(i),y2(i)],'r')
end
axis equal
axis([-1,1,-1,1])
grid
title('x-y-plane')
xlabel('x(m)')
ylabel('y(m)')
% legend('joint_1','end point')
energy_cost
unweighted_terminal_cost
end
function xf = rk4(ode,h,t,x,u)
  k1 = ode(t,x,u);
  k2 = ode(t,x+h/2*k1,u);
  k3 = ode(t,x+h/2*k2,u);
  k4 = ode(t,x+h*k3,  u);
  xf = x + h/6 * (k1 + 2*k2 + 2*k3 + k4); 
end
function xf = euler(ode,h,t,x,u)
  k1 = ode(t,x,u);
  xf = x + h * (k1); 
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

