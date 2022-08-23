function disturbance_compensation_nonlinear_test()
clear
[A,B,B_dbth,C,D]=symbolic_rotor();
A_drone = A;
A_drone([4 8],:)=[];
% A_drone(8,:)=[];
A_drone(:,[4 8])=[];
% A_drone(:,8)=[];
B_drone = B_dbth;
B_drone([4 8],:)=[];
% B_drone(8,:)=[];
B_drone(:,3)=[];
E_drone=[A(:,4),B(:,3)];
E_drone([4 8],:)=[];
C_drone = C;
C_drone(3,:)=[];
C_drone(:,[4 8])=[];
% E(8,:)=[];
K_z_drone = B_drone\E_drone;
% K_z_drone=zeros(2,2);
A_arm=[0,1;A(8,4),0];
B_arm=[0;B_dbth(8,3)];
E_arm=[0;A(8,3)];
K_z_arm = B_arm\E_arm;
% K_z_arm = 0;
C_arm=[1,0];

Q_drone=diag([5,25,0.1, 1,3,1]);
% more for x,z: faster towards desired pos, but with overshoot
% more for x_dot, z_dot; slower with no overshoot
% more for beta: beta less affected by positioning
% more for beta_dot: beta tremble less
% more for theta_dot: less difference of motors, x slower
R_drone=diag([100,100]);
% more f1, f3: slow down
% more tau: worse beta for lowering cost of tau, less peak of tau
%% calculation of the solution of the riccati equation and the controller

[P_drone,L_drone,K_drone] = care(A_drone,B_drone,Q_drone,R_drone,zeros(size(B_drone)),eye(size(A_drone)));
V_drone=inv(C_drone*((B_drone*K_drone-A_drone)\B_drone));

Q_arm=diag([20, 1]);
R_arm=diag([10]);

%% calculation of the solution of the riccati equation and the controller

[P_arm,L_arm,K_arm] = care(A_arm,B_arm,Q_arm,R_arm,zeros(size(B_arm)),eye(size(A_arm)));
V_arm=inv(C_arm*((B_arm*K_arm-A_arm)\B_arm));
% te=10;
% 
% x0=zeros(8,1);
% sim('linearized_anti_disturbance',te)
% figure(2)
% subplot(2,1,1)
% plot(decoupled_output.time, [decoupled_output.signals.values],decoupled_setpoints.time,[decoupled_setpoints.signals.values],'--')
% grid
% legend('x','z','beta')
% xlabel('t(s)')
% ylabel('position (m), angle (rad)')
% title('Decoupled Output')
% subplot(2,1,2)
% plot(decoupled_input.time, [decoupled_input.signals.values])
% grid
% legend('delta f1','delta f3','tau')
% xlabel('t(s)')
% ylabel('force (N), torque (Nm)')
% title('Decoupled Input')
% % C_theta = [0,0,1,0, 0,0,0,0];
% % sim('linearized_model',te)
% % sys = ss(A-B*K,V,C,D)
% % impulse(sys)
N=1000;
tF=10;
t=linspace(tF/N,tF,N);
x0=[0;0;0;pi/2;0;0;0;0];
[mg_num, mq_num, Ig_num, Iq_num, L_center_of_mass, ~, gravity_num] = project_parameters(0);
u0 = [(mg_num+mq_num)*gravity_num/2;(mg_num+mq_num)*gravity_num/2;0];
x_step = 1.0*double(t>=1);
z_step = 0.0*double(t>=1);
beta_step = 0.6*double(t>=6);
delta_x_new = x0-x0;
x_new = x0;
u_new =u0;
delta_u_new = zeros(3,1);
figure(3)
subplot(3,1,1)
hold on
subplot(3,1,2)
hold on
subplot(3,1,3)
hold on
for i=1:N
    x_old = x_new;
    u_old = u_new;
    delta_ud=V_drone*[x_step(i);z_step(i)]-K_drone*delta_x_new([1 2 3 5 6 7])-K_z_drone*[delta_x_new(4);delta_u_new(3)];
    delta_ua = V_arm*beta_step(i)-K_arm*delta_x_new([4 8])-K_z_arm*delta_x_new(3);
    delta_u_new=[delta_ud;delta_ua];
    u_new = delta_u_new+u0;
    u_f_M = [(u_new(1)+u_new(2));-u_new(1)+u_new(2);u_new(3)];
    x_new = rk4(@(t,x,u)diff_eq_grasping(t,x,u),tF/N,0,x_new,u_f_M);
    delta_x_new = x_new-x0;
    subplot(3,1,1)
    plot([t(i)-tF/N;t(i)],[x_old(1);x_new(1)],'b')
    hold on
    plot([t(i)-tF/N;t(i)],[x_old(2);x_new(2)],'g')
    plot([t(i)-tF/N;t(i)],[x_old(4);x_new(4)]-pi/2,'k')
    subplot(3,1,2)
    plot([t(i)-tF/N;t(i)],[u_old(1);u_new(1)],'b')
    hold on
    plot([t(i)-tF/N;t(i)],[u_old(2);u_new(2)],'g')
    subplot(3,1,3)
    plot([t(i)-tF/N;t(i)],[u_old(3);u_new(3)],'b')
    hold on
end
subplot(3,1,1)
plot(t,x_step,'--',t,z_step,'--',t,beta_step,'--')
hold off
legend('x','z','beta')
xlabel('t(s)')
ylabel('position (m), angle (rad)')
title('Decoupled Output')
grid
subplot(3,1,2)
hold off
grid
legend('delta f1','delta f3')
xlabel('t(s)')
ylabel('force (N)')
title('Decoupled Input rotors')
subplot(3,1,3)
hold off
grid
legend('tau')
xlabel('t(s)')
ylabel('torque (Nm)')
title('Decoupled Input tau')
end
function dx=diff_eq_grasping(t,x,u)
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
[mg, mq, Ig, Iq, L, ~, gravity]=project_parameters(0);
D = [mg+mq,0,0,-mg*L*sin(beta);0,mg+mq,0,-mg*L*cos(beta);0,0,Iq,0;-mg*L*sin(beta),-mg*L*cos(beta),0,mg*L*L+Ig];
F=[thrust*sin(theta);thrust*cos(theta);M_drone-M_arm;M_arm];
C=[0,0,0,-L*mg*cos(beta)*beta_dot;0,0,0,L*mg*sin(beta)*beta_dot;0,0,0,0;0,0,0,0];
G=[0;(mg+mq)*gravity;0;-mg*gravity*L*cos(beta)];
q_dot = x(5:8);
q_dotdot=D\(F-C*q_dot-G);
dx=[q_dot;q_dotdot];
end
function xf = rk4(ode,h,t,x,u)
  k1 = ode(t,x,u);
  k2 = ode(t,x+h/2*k1,u);
  k3 = ode(t,x+h/2*k2,u);
  k4 = ode(t,x+h*k3,  u);
  xf = x + h/6 * (k1 + 2*k2 + 2*k3 + k4); 
end