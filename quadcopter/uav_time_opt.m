function uav_time_opt()
addpath('../../../casadi')
import casadi.*;
clear;
clc;

%% Parameters and constraints
[gravity,mass,kF,kM,L] = project_parameters;
x0pos=zeros(3,1);
x0vel=zeros(3,1);
x0angle=[0;0;3];
x1pos=[2;0;1.5];
x1vel=[nan;0;0];
x2pos=[6;4;1.5];
x2vel=[0;nan;0];
xFpos=[7;6;0.2];
xFvel=zeros(3,1);
xFangle=[0;0;3];

figure(3)
hold on
xlabel('x')
ylabel('y')
zlabel('z')
axis equal
grid on
x1vel_normalized=[1;0;0];
x2vel_normalized=[0;1;0];
plot_rings(0.3,atan2(x1vel_normalized(2),x1vel_normalized(1)),x1pos(1),x1pos(2),x1pos(3),3)
plot_rings(0.3,atan2(x2vel_normalized(2),x2vel_normalized(1)),x2pos(1),x2pos(2),x2pos(3),3)
% R=diag([1e-6,1e-6,1e-6,1e-6]);
R=eye(4);
vel_constr = 5;
angle_constr = pi/4;
input_constr = 6;

%% There are two way points, the trajectory is divided into 3 interval pieces,
%% and then discreticized into 301 sections, 100 for each piece
N=100;
opti=casadi.Opti();
X=opti.variable(12,3*N+1);
U=opti.variable(4,3*N);
t_interval=opti.variable(3,1);

%% Setting the objective and constraints for the first interval
opti.subject_to(X(1:3,1) == x0pos);
opti.subject_to(X(4:6,1) == x0vel);
opti.subject_to(X(7:9,1) == x0angle);
tau=0;
J1=0;
for i=1:N
    tau=tau+1/N;
    x_next = rk4(@(t,x,u)diff_eq_time_trans(t,x,u,t_interval(1)),1/N,tau,X(:,i),U(:,i));
    opti.subject_to(X(:,i+1)==x_next);
    J1 = J1 + U(:,i)' * R * U(:,i)*1/N;
end
opti.subject_to(X(10:12,1) == zeros(3,1));
opti.subject_to(X(1:3,N+1) == x1pos);
opti.subject_to(X(5:6,N+1) == x1vel(2:3));
%% Setting the objective and constraints for the second interval
J2=0;
for i=N+1:2*N
    tau=tau+1/N;
    x_next = rk4(@(t,x,u)diff_eq_time_trans(t,x,u,t_interval(2)),1/N,tau,X(:,i),U(:,i));
    opti.subject_to(X(:,i+1)==x_next);
    J2 = J2 + U(:,i)' * R * U(:,i)*1/N;
end
opti.subject_to(X(1:3,2*N+1) == x2pos);
opti.subject_to(X(4,2*N+1) == x2vel(1));
opti.subject_to(X(6,2*N+1) == x2vel(3));
%% Setting the objective and constraints for the third interval
J3=0;
for i=2*N+1:3*N
    tau=tau+1/N;
    x_next = rk4(@(t,x,u)diff_eq_time_trans(t,x,u,t_interval(3)),1/N,tau,X(:,i),U(:,i));
    opti.subject_to(X(:,i+1)==x_next);
    J3 = J3 + U(:,i)' * R * U(:,i)*1/N;
end
opti.subject_to(X(1:3,end) == xFpos);
opti.subject_to(X(4:6,end) == xFvel);
opti.subject_to(X(7:9,end) == xFangle);
opti.subject_to(X(10:12,end) == zeros(3,1));
%% Bounds for states and input
opti.subject_to(-vel_constr <= X(4:6,:) <= vel_constr);
opti.subject_to(-angle_constr <= X(7:8,:) <= angle_constr);
opti.subject_to(0 <= U(1:4,:) <= input_constr);
opti.subject_to(X(9,:) == 3);
opti.subject_to(t_interval(1)>= 0);
opti.subject_to(t_interval(2)>= 0);
opti.subject_to(t_interval(3)>= 0);
%% add up the objective
% opti.minimize(t_interval(1)+t_interval(2)+t_interval(3)+J1+J2+J3);
% opti.minimize(J1+J2+J3);
opti.minimize(J1*t_interval(1)+J2*t_interval(2)+J3*t_interval(3));
%% choose the solver
opti.solver('ipopt');
%% Give initial guess, without this, it will take far longer time to solve, or
%% sometimes the solution cannot be found. 
%% Initialize the input values to be the the hovering thrusts
opti.set_initial(t_interval,[1;1;1]);
opti.set_initial(U,0.25*gravity*mass*ones(4,3*N));

sol = opti.solve();
%% collect the solved data
pos=sol.value(X(1:3,:))';
vel=sol.value(X(4:6,:))';
angle=sol.value(X(7:9,:))';
angvel=sol.value(X(10:12,:))';
rotor=sol.value(U(1:4,:))';
%% organize the time axes
tx1=linspace(0,sol.value(t_interval(1)),N+1);
tx2=linspace(0,sol.value(t_interval(2)),N+1)+tx1(end);
tx3=linspace(0,sol.value(t_interval(3)),N+1)+tx2(end);
tx=[tx1,tx2(2:end),tx3(2:end)];
tu=tx(1:end-1);
%% plotting
figure(1)
subplot(4,1,1)
plot(tx,pos)
hold on
plot(tx(length(tx1)),pos(length(tx1)),'*')
plot(tx(length(tx)-length(tx3)),pos(length(tx)-length(tx3)),'*')
hold off
% title('pos')
xlabel('t (s)')
ylabel('position (m)')
legend('x','y','z')
grid
subplot(4,1,2)
plot(tx,vel)
% title('vel')
xlabel('t (s)')
ylabel('velocity (m/s)')
legend('x','y','z')
grid
subplot(4,1,3)
plot(tx,angle)
% title('acc')
xlabel('t (s)')
ylabel('angle (rad)')
legend('x','y','z')
grid
subplot(4,1,4)
plot(tu,rotor)
% title('rotors')
xlabel('t (s)')
ylabel('rotor Forces(N)')
legend('1','2','3','4')
grid
figure(2)
plot(tx,angvel)
xlabel('t (s)')
ylabel('angular velocity (rad/s)')
legend('x','y','z')
grid
figure(3)
plot3(pos(:,1),pos(:,2),pos(:,3))


end
function xf = rk4(ode,h,t,x,u)
  k1 = ode(t,x,u);
  k2 = ode(t,x+h/2*k1,u);
  k3 = ode(t,x+h/2*k2,u);
  k4 = ode(t,x+h*k3,  u);
  xf = x + h/6 * (k1 + 2*k2 + 2*k3 + k4); 
end
function dx=diff_eq_time_trans(tau,x,u,tf)
[gravity,mass,kF,kM,L] = project_parameters;
pos=x(1:3);
vel=x(4:6);
angle=x(7:9);
dangle_dt=x(10:12);
dangvel_B_dt=[kF*L*(u(2)-u(4));kF*L*(-u(1)+u(3));kM*(u(1)-u(2)+u(3)-u(4))];
abs_acc = (u(1)+u(2)+u(3)+u(4))/mass;
cp = cos(angle(2));
sp = sin(angle(2));
sr = sin(angle(1));
cr = cos(angle(1));
sy = sin(angle(3));
cy = cos(angle(3));
% % 	R_B2W(1,1) = cp * cy;
% % 	R_B2W(1,2) = ((sr * sp * cy) - (cr * sy));
% R_B2W_13 = ((cr * sp * cy) + (sr * sy));
% % 	R_B2W(2,1) = cp * sy;
% % 	R_B2W(2,2) = ((sr * sp * sy) + (cr * cy));
% R_B2W_23 = ((cr * sp * sy) - (sr * cy));
% % 	R_B2W(3,1) = -sp;
% % 	R_B2W(3,2) = sr * cp;
% R_B2W_33 = cr * cp;	
  R_B2W_11 = cp * cy;
  R_B2W_21 = ((sr * sp * cy) - (cr * sy));
  R_B2W_31 = ((cr * sp * cy) + (sr * sy));
  R_B2W_12 = cp * sy;
  R_B2W_22 = ((sr * sp * sy) + (cr * cy));
  R_B2W_32 = ((cr * sp * sy) - (sr * cy));
  R_B2W_13 = -sp;
  R_B2W_23 = sr * cp;
  R_B2W_33 = cr * cp; 
new_acc = abs_acc*[R_B2W_13;R_B2W_23;R_B2W_33]-[0;0;gravity];
angrate_W_x = R_B2W_11*dangle_dt(1) + R_B2W_12*dangle_dt(2) + R_B2W_13*dangle_dt(3);
angrate_W_y = R_B2W_21*dangle_dt(1) + R_B2W_22*dangle_dt(2) + R_B2W_23*dangle_dt(3);
angrate_W_z = R_B2W_31*dangle_dt(1) + R_B2W_32*dangle_dt(2) + R_B2W_33*dangle_dt(3);
dangle_dt=[angrate_W_x;angrate_W_y;angrate_W_z];
dvel_dt = new_acc;
dpos_dt = vel;
dx=tf*[dpos_dt;dvel_dt;dangle_dt;dangvel_B_dt];
end
function dx=diff_eq(t,x,u)
[gravity,mass,kF,kM,L] = project_parameters;
pos=x(1:3);
vel=x(4:6);
angle=x(7:9);
% yaw=x(10);
dangle_dt=x(10:12);

% zB=1/(norm(acc))*acc;
% xC=[cos(yaw);sin(yaw);0];
% yB=cross(zB,xC);
% yB=1/(norm(yB))*yB;
% xB=cross(yB,zB);
% R_B2W = [xB,yB,zB];
% quat_B2W=rotm2quat(R_B2W);
% q1=quat_B2W(1);
% q2=quat_B2W(2);
% q3=quat_B2W(3);
% q4=quat_B2W(4);

dangvel_B_dt=[kF*L*(u(2)-u(4));kF*L*(-u(1)+u(3));kM*(u(1)-u(2)+u(3)-u(4))];
% quat_update_matrix = 0.5*[0,-angvel_B(1),-angvel_B(2),-angvel_B(3);
%     angvel_B(1),0,angvel_B(3),-angvel_B(2);
%     angvel_B(2),-angvel_B(3),0,angvel_B(1);
%     angvel_B(3),angvel_B(2),angvel_B(1),0];
% q_dot=quat_update_matrix*quat_B2W;
% q_dot1=q_dot(1);
% q_dot2=q_dot(2);
% q_dot3=q_dot(3);
% q_dot4=q_dot(4);
abs_acc = (u(1)+u(2)+u(3)+u(4))/mass;
R_B2W=eye(3);
	 cp = cos(angle(2));
	 sp = sin(angle(2));
	 sr = sin(angle(1));
	 cr = cos(angle(1));
	 sy = sin(angle(3));
	 cy = cos(angle(3));
	R_B2W(0,0) = cp * cy;
	R_B2W(0,1) = ((sr * sp * cy) - (cr * sy));
	R_B2W(0,2) = ((cr * sp * cy) + (sr * sy));
	R_B2W(1,0) = cp * sy;
	R_B2W(1,1) = ((sr * sp * sy) + (cr * cy));
	R_B2W(1,2) = ((cr * sp * sy) - (sr * cy));
	R_B2W(2,0) = -sp;
	R_B2W(2,1) = sr * cp;
	R_B2W(2,2) = cr * cp;	
new_acc = abs_acc*R_B2W(:,3)-[0;0;gravity];
% dzB_dt=[2*q_dot2*q4+2*q2*q_dot4+2*q_dot1*q3+2*q1*q_dot3;
%     2*q_dot3*q4+2*q3*q_dot4-2*q_dot1*q2-2*q1*q_dot2;
%     -4*q_dot2*q2-4*q_dot3*q3];
% dacc_dt = abs_acc * dzB_dt;
dvel_dt = new_acc;
dpos_dt = vel;
% dyaw_dt = angvel_B(3);
dx=[dpos_dt;dvel_dt;dangle_dt;dangvel_B_dt];
end
function [gravity,mass,kF,kM,L] = project_parameters
%% Definition of system parameters
kF=10;
kM=2;
mass=1.0;
L=0.225;
gravity=9.81;
end

function plot_rings(radius,yaw,x_tr,y_tr,z_tr,figure_num)
% Original points, original plane
t = linspace(0,2*pi);
x = radius*cos(t);
y = radius*sin(t);
z = 0*t;
pnts = [x;y;z];
% unit normal for original plane
n0 = [0;0;1]; 
n0 = n0/norm(n0);
% unit normal for plane to rotate into 
% plane is orthogonal to n1... given by equation
% n1(1)*x + n1(2)*y + n1(3)*z = 0
n1 = [cos(yaw);sin(yaw);0]; 
n1 = n1/norm(n1); 
% theta is the angle between normals
c = dot(n0,n1) / ( norm(n0)*norm(n1) ); % cos(theta)
s = sqrt(1-c*c);                        % sin(theta)
u = cross(n0,n1) / ( norm(n0)*norm(n1) ); % rotation axis...
u = u/norm(u); % ... as unit vector
C = 1-c;
% the rotation matrix
R = [u(1)^2*C+c, u(1)*u(2)*C-u(3)*s, u(1)*u(3)*C+u(2)*s
    u(2)*u(1)*C+u(3)*s, u(2)^2*C+c, u(2)*u(3)*C-u(1)*s
    u(3)*u(1)*C-u(2)*s, u(3)*u(2)*C+u(1)*s, u(3)^2*C+c];
% Rotated points
newPnts = R*pnts;
newPnts(1,:)=newPnts(1,:)+x_tr;
newPnts(2,:)=newPnts(2,:)+y_tr;
newPnts(3,:)=newPnts(3,:)+z_tr;
plot3(newPnts(1,:),newPnts(2,:),newPnts(3,:),'ko') 
hold on
z_stem=0:0.05:(z_tr-radius);
x_stem=ones(1,length(z_stem))*x_tr;
y_stem=ones(1,length(z_stem))*y_tr;
figure(figure_num)
plot3(x_stem,y_stem,z_stem,'ko')
axis equal
end