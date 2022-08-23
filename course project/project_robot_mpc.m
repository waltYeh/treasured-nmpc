function project_robot_mpc()
addpath('../../../../casadi')
import casadi.*;
clear;
clc;
control_interval=0.2;
noise_var=0.1;
x0=[-5;-4;0;0];
xF=[pi/2;0;0;0];
tF=3;
t=0;
x=x0;
figure(1)
control_loop_cnt=0;
while 1
    control_loop_cnt=control_loop_cnt+1;
    est=estimation(x,noise_var);
    t_est=t;
    [u,tu]=ocp(t, est, tF, xF);
    t_span=[t, t+control_interval];
    
    [tode,xode] = ode45(@(t,x)diff_eq_for_ode(t,x,u,tu),t_span,x);
    t=tode(end);
    x=xode(end,:)';

    subplot(4,1,1)
    plot(tode,xode(:,1),'b',tode,xode(:,2),'r')
    axis([0 3.1 -6 5])
    hold on
    plot(t_est,est(1),'b*')
    plot(t_est,est(2),'r*')
    title('states')
    xlabel('t (s)')
    ylabel('angle (rad)')
    legend('q_1','q_2')

    subplot(4,1,2)
    plot(tode,xode(:,3),'b',tode,xode(:,4),'r')
    axis([0 3.1 -5 5])
    hold on
    plot(t_est,est(3),'b*')
    plot(t_est,est(4),'r*')
    title('states')
    xlabel('t (s)')
    ylabel('angular velocity (rad/s)')
    legend('omiga_1','omiga_2')

    subplot(4,1,3)
%     tustairs=[tu',tu'];
%     ustairs=[u(1,:)',u(2,:)'];
    stairs(tu',u(1,:)','b')
    hold on
    stairs(tu',u(2,:)','r')
%     plot(tu,u(1,:),'b',tu,u(2,:),'r')
    axis([0 3.1 -1010 1010])
    title('optimized inputs for each step')
    xlabel('t (s)')
    ylabel('torque (Nm)')
    legend('u_1','u_2')

    subplot(4,1,4)
    index=max(find(tu<=t+0.001));
    stairs(tu(1:index)',u(1,1:index)','b')
    hold on
    stairs(tu(1:index)',u(2,1:index)','r')
%     plot(tu(1:index),u(1,1:index),'b',tu(1:index),u(2,1:index),'r')
    axis([0 3.1 -1010 1010])
    
    title('real inputs into the system')
    xlabel('t (s)')
    ylabel('torque (Nm)')
    legend('u_1','u_2')
    if t>=tF%-control_interval
        break;
    end
    
end
error=x-xF
end
%%
function [u,tu]=ocp(tk, est0, tf, xf)
persistent run_cnt
if isempty(run_cnt)
    run_cnt=0;
end
run_cnt=run_cnt+1;
opti=casadi.Opti();

opt_discrete=0.025;
Q=diag([0,0,0,0]);
R=diag([10,10]);
W=5e8;%*(tf-tk)/3;
if tf-tk<0.5
    N=round((1)/opt_discrete);
else
    N=round((tf-tk)/opt_discrete);
end
[~, ~, ~, ~, state_constr, input_constr] = project_parameters;
X=opti.variable(4,N+1);
U=opti.variable(2,N);

J=0;
for i=1:N
    x_next = rk4(@(t,x,u)diff_eq_for_ocp(t,x,u),opt_discrete,0,X(:,i),U(:,i));
    opti.subject_to(X(:,i+1)==x_next);
    error = x_next - xf;
    J = J + (U(:,i)' * R * U(:,i)+error'*Q*error)*opt_discrete;
end
opti.subject_to(X(:,1) == est0);
% if W==0
%     opti.subject_to(X(:,end) == xf);
% end
J = J + (X(:,end) - xf)'*(X(:,end) - xf)*W;%*(tf-tk)/tf;
opti.subject_to(-state_constr <= X(3,:) <= state_constr);
opti.subject_to(-state_constr <= X(4,:) <= state_constr);
opti.subject_to(-input_constr <= U(1,:) <= input_constr);
opti.subject_to(-input_constr <= U(2,:) <= input_constr);
if tk<tf/4
%     i=round(tf/6/opt_discrete);
%     opti.subject_to(X(3,i)>=X(4,i));
%     i=N-round(tf/6/opt_discrete);
%     opti.subject_to(X(3,i)<=X(4,i));
    free_zone=round(tf/12/opt_discrete);
    i=round(tf/6/opt_discrete);
    opti.subject_to(X(3,free_zone:i)>=X(4,free_zone:i));
    i=N-round(tf/6/opt_discrete);
    opti.subject_to(X(3,i:end-free_zone)<=X(4,i:end-free_zone));

end
opti.minimize(J);
opti.solver('ipopt');
sol = opti.solve();
u1=sol.value(U(1,:));
u2=sol.value(U(2,:));
tu=linspace(tk,tf-opt_discrete,N);
u=[u1;u2];

%     subplot(4,1,1)
%     plot([tu,tf],sol.value(X(1,:)),':.g',[tu,tf],sol.value(X(2,:)),':.g')
%     
%     hold on
%     
% 
%     subplot(4,1,2)
%     plot([tu,tf],sol.value(X(3,:)),':.g',[tu,tf],sol.value(X(4,:)),':.g')
%     hold on
end
%%
function est=estimation(x,noise_var)

est_q1=x(1)+noise_var*randn;
est_q2=x(2)+noise_var*randn;
est_q1_dot=x(3)+noise_var*randn;
est_q2_dot=x(4)+noise_var*randn;
est=[est_q1;est_q2;est_q1_dot;est_q2_dot];
end
function xf = rk4(ode,h,t,x,u)
  k1 = ode(t,x,u);
  k2 = ode(t,x+h/2*k1,u);
  k3 = ode(t,x+h/2*k2,u);
  k4 = ode(t,x+h*k3,  u);
  xf = x + h/6 * (k1 + 2*k2 + 2*k3 + k4); 
end
%%
function dx=diff_eq_for_ocp(t,x,u)
q1=x(1);
q2=x(2);
w1=x(3);
w2=x(4);
[b,c,g,~,~,~]=project_parameters;
% b(1)=180;
% b(2)=45;
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
%%
function dx=diff_eq_for_ode(t,x,U,tU)
u1=interp1(tU,U(1,:),t,'previous','extrap');
u2=interp1(tU,U(2,:),t,'previous','extrap');
u=[u1;u2];
if isnan(u(1))
    keyboard
end
if isnan(u(2))
    keyboard
end
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

% cond(B)
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

