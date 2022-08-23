function uav_grasp_opt()
addpath('../../../casadi')
import casadi.*;

clear;
clc;
x0=[-3;1;0;0;0;0];
xF=[6;1;0;0;0;0];
x_object = 1;
z_object = 1;
R=diag([1;10]);
vel_constr = 5;
angle_constr = pi/10;
thrust_constr = 20;
M_drone_constr = 5*0.34/2;
[mg, mq, Ig, Iq, L_center_of_mass, L_arm, gravity] = project_parameters;
%% There are two way points, the trajectory is divided into 3 interval pieces,
%% and then discreticized into 301 sections, 100 for each piece
N=100;
opti=casadi.Opti();
X=opti.variable(6,2*N+1);
U=opti.variable(2,2*N);
t_interval=opti.variable(2,1);

%% Setting the objective and constraints for the first interval
opti.subject_to(X(:,1) == x0);
tau=0;
J1=0;
for i=1:N
    tau=tau+1/N;
    x_next = rk4(@(t,x,u)diff_eq_grasping_time_trans(t,x,u,t_interval(1)),1/N,tau,X(:,i),U(:,i));
    opti.subject_to(X(:,i+1)==x_next);
    J1 = J1 + U(:,i)' * R * U(:,i)*1/N;
end
% end point of arm overlaps with the object for grasping, vel should be 0
opti.subject_to(X(1,N+1) == x_object);
opti.subject_to(X(2,N+1) == z_object);
opti.subject_to(X(4,N+1) == 0);
opti.subject_to(X(5,N+1) == 0);
%% Setting the objective and constraints for the second interval
J2=0;
for i=N+1:2*N
    tau=tau+1/N;
    x_next = rk4(@(t,x,u)diff_eq_grasping_time_trans(t,x,u,t_interval(2)),1/N,tau,X(:,i),U(:,i));
    opti.subject_to(X(:,i+1)==x_next);
    J2 = J2 + U(:,i)' * R * U(:,i)*1/N;
end
opti.subject_to(X(:,end) == xF);
%% Bounds for states and input
opti.subject_to(-vel_constr <= X(4:5,:) <= vel_constr);
opti.subject_to(-angle_constr <= X(3,:) <= angle_constr);
opti.subject_to(0 <= U(1,:) <= thrust_constr);
opti.subject_to(-M_drone_constr <= U(2,:) <= M_drone_constr);
opti.subject_to(X(2,:) == z_object);
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
U_init = zeros(2,2*N);
U_init(1,:) = 0.5*gravity*(mg+mq);
opti.set_initial(U,U_init);

sol = opti.solve();
%% collect the solved data
x_opt=sol.value(X(1,:))';
z_opt=sol.value(X(2,:))';
theta_opt=sol.value(X(3,:))';

xdot_opt=sol.value(X(4,:))';
zdot_opt=sol.value(X(5,:))';
thetadot_opt=sol.value(X(6,:))';

thrust_opt=sol.value(U(1,:))';
M_drone_opt=sol.value(U(2,:))';
[xg_opt,zg_opt] = state_arm_end(x_opt,z_opt,theta_opt,L_arm);

%% organize the time axes
tx1=linspace(0,sol.value(t_interval(1)),N+1);
tx2=linspace(0,sol.value(t_interval(2)),N+1)+tx1(end);
tx=[tx1,tx2(2:end)];
tu=tx(1:end-1);
%% plotting
figure(1)
subplot(4,1,1)
plot(tx,x_opt,tx,z_opt)
% title('pos')
xlabel('t (s)')
ylabel('pos (m)')
legend('x','z')
grid
subplot(4,1,2)
plot(tx, xdot_opt,tx,zdot_opt)
% title('vel')
xlabel('t (s)')
ylabel('vel (m/s)')
legend('x','z')
grid
subplot(4,1,3)
plot(tx,theta_opt)
% title('acc')
xlabel('t (s)')
ylabel('angle (rad)')
legend('theta')
grid
subplot(4,1,4)
plot(tx,thetadot_opt)
% title('acc')
xlabel('t (s)')
ylabel('anglar rate (rad/s)')
legend('theta')
figure(2)
subplot(2,1,1)
plot(tu,thrust_opt)
xlabel('t (s)')
ylabel('thrust (N)')
% legend('x','y','z')
grid
subplot(2,1,2)
plot(tu,M_drone_opt)
xlabel('t (s)')
ylabel('M drone (Nm)')
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
    if mod(i,5)==1
        plot([x_opt(i),xg_opt(i)],[z_opt(i),zg_opt(i)],'-o')
        quiver(x_opt(i),z_opt(i),cos(theta_opt(i)),-sin(theta_opt(i)),0.05)
    end
end
hold off
grid
dx_dt = zeros(6,2*N+1);
X_all_opt = sol.value(X);
U_all_opt = sol.value(U);
for i=2:length(tu)
    dx_dt(:,i) = diff_eq_grasping_time_trans(0,X_all_opt(:,i),U_all_opt(:,i),1);
end
figure(4)
plot(tx,dx_dt(4,:))
xlabel('t (s)')
ylabel('acc (m/s^2)')
legend('x')
grid

acc_opt = dx_dt(4,:)';
M_opt = [M_drone_opt(1);M_drone_opt];
Omg_opt = X_all_opt(6,:)';
assignin('base','acc',acc_opt);
assignin('base','M',M_opt);
assignin('base','w',Omg_opt);
assignin('base','N',N);
Ts = tx(end)/2/N;
assignin('base','T',Ts);
I=mg*L_center_of_mass*L_center_of_mass+Ig+Iq;
assignin('base','I',I);
assignin('base','mg',mg);
assignin('base','gravity',gravity);
assignin('base','L',L_center_of_mass);

% y=zeros(2*N-2,1);
% M=zeros(2*N-2,3);
% for k = 4:2*N+1
%     y(k-2)=Omg_opt(k)+Omg_opt(k-2);
%     M(k-2,1) = M_opt(k)-M_opt(k-2);
%     M(k-2,2) = acc_opt(k)-acc_opt(k-2);
%     M(k-2,3) = -Omg_opt(k-1);
% end
% theta0=M\y;
end
function [xg,zg] = state_arm_end(x,z,theta,L_arm)
xg=x+L_arm*cos(theta+pi/2);
zg=z-L_arm*sin(theta+pi/2);
end
function dx=diff_eq_grasping_time_trans(t,x,u,t_trans)
% xq=x(1);
% zq=x(2);
theta=x(3);
% xq_dot=x(5);
% zq_dot=x(6);
theta_dot=x(6);
thrust = u(1);
M_drone = u(2);
[mg, mq, Ig, Iq, L, ~, gravity]=project_parameters;
D = [mg+mq,0,-mg*L*cos(theta);0,mg+mq,mg*L*sin(theta);-mg*L*cos(theta),mg*L*sin(theta),mg*L*L+Ig+Iq];
F=[thrust*sin(theta);thrust*cos(theta);M_drone];
C=[0,0,L*mg*sin(theta)*theta_dot;0,0,L*mg*cos(theta)*theta_dot;0,0,0];
G=[0;(mg+mq)*gravity;-mg*gravity*L*sin(theta)];
q_dot = x(4:6);
q_dotdot=D\(F-C*q_dot-G);
dx=[q_dot;q_dotdot]*t_trans;
end
function [mg, mq, Ig, Iq, L_center_of_mass, L_arm, gravity] = project_parameters
%% Definition of system parameters
mq = 0.531;
mg = 0.2;
h=0.02;
L_center_of_mass=0.1;
Ig = 1/12*mg*(h^2+L_center_of_mass^2);
Iq = 0.00368;
% http://wiki.asctec.de/display/AR/CAD+Models
gravity = 9.81;
L_arm = 0.2;
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