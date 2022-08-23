function uav_grasp_opt_var_mass()
addpath('../../../casadi')
import casadi.*;

clear;
clc;
payload = 0.12;
x0=[-3;1.5;0;pi/2;0;0;0;0];
xF=[3;1;0;pi/2;0;0;0;0];
x_object = 0;
z_object = 0;
R=diag([1;1;10]);
Q_beta = 5;
vel_constr = 5;
angle_constr = pi/4;
arm_range_min = 0;
arm_range_max = pi;
thrust_constr = 20;
M_drone_constr = 5*0.34/2;
M_arm_constr = 3;
[mg, mq, ~, ~, ~, L_arm, gravity] = project_parameters(0);
%% The trajectory is divided into 2 interval pieces, before grasping and after
%% and then discreticized into 201 sections, 100 for each piece
N=100;
opti=casadi.Opti();
X=opti.variable(8,2*N+1);
U=opti.variable(3,2*N);
t_interval=opti.variable(2,1);

%% Setting the objective and constraints for the first interval
opti.subject_to(X(:,1) == x0);
tau=0;
J1=0;
for i=1:N
    tau=tau+1/N;
    x_next = rk4(@(t,x,u)diff_eq_grasping_time_trans(t,x,u,t_interval(1),0),1/N,tau,X(:,i),U(:,i));
    opti.subject_to(X(:,i+1)==x_next);
    J1 = J1 + U(:,i)' * R * U(:,i)*1/N + (X(4,i)-X(3,i)-pi/2)'* Q_beta * (X(4,i)-X(3,i)-pi/2)*1/N;
end
% end point of arm overlaps with the object for grasping, vel should be 0
opti.subject_to(X(1,N+1)+L_arm*cos(X(4,N+1)) == x_object);
opti.subject_to(X(2,N+1)-L_arm*sin(X(4,N+1)) == z_object);
opti.subject_to(X(5,N+1)-L_arm*sin(X(4,N+1))*X(8,N+1) == 0);
opti.subject_to(X(6,N+1)-L_arm*cos(X(4,N+1))*X(8,N+1) == 0);
%% Setting the objective and constraints for the second interval
J2=0;
for i=N+1:2*N
    tau=tau+1/N;
    x_next = rk4(@(t,x,u)diff_eq_grasping_time_trans(t,x,u,t_interval(2),payload),1/N,tau,X(:,i),U(:,i));
    opti.subject_to(X(:,i+1)==x_next);
    J2 = J2 + U(:,i)' * R * U(:,i)*1/N + (X(4,i)-X(3,i)-pi/2)'* Q_beta * (X(4,i)-X(3,i)-pi/2)*1/N;
end
opti.subject_to(X(:,end) == xF);
%% Bounds for states and input
opti.subject_to(-vel_constr <= X(5:6,:) <= vel_constr);
opti.subject_to(-angle_constr <= X(3,:) <= angle_constr);
opti.subject_to(arm_range_min <= X(4,:)-X(3,:) <= arm_range_max);
opti.subject_to(0 <= U(1,:) <= thrust_constr);
opti.subject_to(-M_drone_constr <= U(2,:) <= M_drone_constr);
opti.subject_to(-M_arm_constr <= U(3,:) <= M_arm_constr);
opti.subject_to(t_interval(1)>= 0);
opti.subject_to(t_interval(2)>= 0);
%% add up the objective
% opti.minimize(t_interval(1)+t_interval(2)+t_interval(3)+J1+J2+J3);
% opti.minimize(J1+J2+J3);
opti.minimize(J1*t_interval(1)+J2*t_interval(2));
%% choose the solver
opti.solver('ipopt');
%% Give initial guess, without this, it will take far longer time to solve, or
%% sometimes the solution cannot be found. 
%% Initialize the input values to be the the hovering thrusts
opti.set_initial(t_interval,[4;4]);
U_init = zeros(3,2*N);
U_init(1,:) = 0.5*gravity*(mg+mq);
opti.set_initial(U,U_init);

sol = opti.solve();
%% collect the solved data
x_opt=sol.value(X(1,:))';
z_opt=sol.value(X(2,:))';
theta_opt=sol.value(X(3,:))';
beta_opt=sol.value(X(4,:))';
xdot_opt=sol.value(X(5,:))';
zdot_opt=sol.value(X(6,:))';
thetadot_opt=sol.value(X(7,:))';
betadot_opt=sol.value(X(8,:))';
thrust_opt=sol.value(U(1,:))';
M_drone_opt=sol.value(U(2,:))';
M_arm_opt=sol.value(U(3,:))';
[xg_opt,zg_opt] = state_arm_end(x_opt,z_opt,beta_opt,L_arm);
%% organize the time axes
tx1=linspace(0,sol.value(t_interval(1)),N+1);
tx2=linspace(0,sol.value(t_interval(2)),N+1)+tx1(end);
tx=[tx1,tx2(2:end)];
tu=tx(1:end-1);
%% plotting
figure(1)
subplot(4,1,1)
plot(tx,x_opt,tx,z_opt)
hold on
plot([tx(N+1),tx(N+1)],[-5,5],'--')
hold off
% title('pos')
xlabel('t (s)')
ylabel('pos (m)')
legend('x','z')
grid
subplot(4,1,2)
plot(tx, xdot_opt,tx,zdot_opt)
hold on
plot([tx(N+1),tx(N+1)],[-vel_constr ,vel_constr ],'--')
hold off
% title('vel')
xlabel('t (s)')
ylabel('vel (m/s)')
legend('x','z')
grid
subplot(4,1,3)
plot(tx,theta_opt,tx,beta_opt)
hold on
plot([tx(N+1),tx(N+1)],[-2,4],'--')
hold off
% title('acc')
xlabel('t (s)')
ylabel('angle (rad)')
legend('theta','beta')
grid
subplot(4,1,4)
plot(tx,thetadot_opt,tx,betadot_opt)
hold on
plot([tx(N+1),tx(N+1)],[-20,40],'--')
hold off
% title('acc')
xlabel('t (s)')
ylabel('anglar rate (rad/s)')
legend('theta','beta')
grid
figure(2)
subplot(3,1,1)
plot(tu,thrust_opt)
hold on
plot([tx(N),tx(N)],[0,thrust_constr],'--')
hold off
xlabel('t (s)')
ylabel('thrust (N)')
% legend('x','y','z')
grid
subplot(3,1,2)
plot(tu,M_drone_opt)
hold on
plot([tx(N),tx(N)],[-M_drone_constr,M_drone_constr],'--')
hold off
xlabel('t (s)')
ylabel('M drone (Nm)')
% legend('x','y','z')
grid
subplot(3,1,3)
plot(tu,M_arm_opt)
hold on
plot([tx(N),tx(N)],[-M_arm_constr,M_arm_constr],'--')
hold off
xlabel('t (s)')
ylabel('M arm (Nm)')
% legend('x','y','z')
grid
figure(3)
plot(x_opt,z_opt)
hold on
plot(xg_opt,zg_opt)
xlabel('x')
ylabel('z')
axis equal
for i=1:length(tx)
    if mod(i,3)==1
        plot([x_opt(i),xg_opt(i)],[z_opt(i),zg_opt(i)],'-o')
        quiver(x_opt(i),z_opt(i),cos(theta_opt(i)),-sin(theta_opt(i)),0.1)
    end
end
hold off
grid
dx_dt = zeros(8,2*N);
X_all_opt = sol.value(X);
U_all_opt = sol.value(U);
for i=1:length(tu)
    if i<=N
        dx_dt(:,i) = diff_eq_grasping_time_trans(0,X_all_opt(:,i),U_all_opt(:,i),1,0);
    else
        dx_dt(:,i) = diff_eq_grasping_time_trans(0,X_all_opt(:,i),U_all_opt(:,i),1,payload);
    end
end
figure(4)
subplot(2,1,1)
plot(tu,dx_dt(5,:),tu,dx_dt(6,:))
hold on
plot([tx(N),tx(N)],[-30,20],'--')
hold off
xlabel('t (s)')
ylabel('acc (m/s^2)')
legend('x','z')
grid
subplot(2,1,2)
plot(tu,dx_dt(7,:),tu,dx_dt(8,:))
hold on
plot([tx(N),tx(N)],[-500,1000],'--')
hold off
xlabel('t (s)')
ylabel('angular acc (rad/s^2)')
legend('theta','beta')
grid
end
function [xg,zg] = state_arm_end(x,z,beta,L_arm)
xg=x+L_arm*cos(beta);
zg=z-L_arm*sin(beta);
end
function dx=diff_eq_grasping_time_trans(t,x,u,t_trans,payload)
% xq=x(1);
% zq=x(2);
theta=x(3);
beta=x(4);
% xq_dot=x(5);
% zq_dot=x(6);
% theta_dot=x(7);
beta_dot=x(8);
thrust = u(1);
M_drone = u(2);
M_arm = u(3);
[mg, mq, Ig, Iq, L, ~, gravity]=project_parameters(payload);
D = [mg+mq,0,0,-mg*L*sin(beta);0,mg+mq,0,-mg*L*cos(beta);0,0,Iq,0;-mg*L*sin(beta),-mg*L*cos(beta),0,mg*L*L+Ig];
F=[thrust*sin(theta);thrust*cos(theta);M_drone-M_arm;M_arm];
C=[0,0,0,-L*mg*cos(beta)*beta_dot;0,0,0,L*mg*sin(beta)*beta_dot;0,0,0,0;0,0,0,0];
G=[0;(mg+mq)*gravity;0;-mg*gravity*L*cos(beta)];
q_dot = x(5:8);
q_dotdot=D\(F-C*q_dot-G);
dx=[q_dot;q_dotdot]*t_trans;
end
function [mg, mq, Ig, Iq, L_center_of_mass, L_arm, gravity] = project_parameters(payload)
%% Definition of system parameters
L_arm = 0.2;
mq = 0.531;
h=0.02;
m_arm = 0.1;
L_center_of_mass=(L_arm/2*m_arm+L_arm*payload)/(m_arm+payload);
mg = m_arm+payload;
Ig = 1/12*mg*(h^2+L_arm^2)+payload*(L_arm/2)^2;
Iq = 0.00368;
% http://wiki.asctec.de/display/AR/CAD+Models
gravity = 9.81;

%% Constraints


%% Scaling factors

end
function xf = rk4(ode,h,t,x,u)
  k1 = ode(t,x,u);
  k2 = ode(t,x+h/2*k1,u);
  k3 = ode(t,x+h/2*k2,u);
  k4 = ode(t,x+h*k3,  u);
  xf = x + h/6 * (k1 + 2*k2 + 2*k3 + k4); 
end