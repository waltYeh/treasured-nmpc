
clear
[A,B,B_dbth,C,D]=symbolic_rotor()
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
te=10;

x0=zeros(8,1);
sim('linearized_anti_disturbance',te)
figure(2)
subplot(2,1,1)
plot(decoupled_output.time, [decoupled_output.signals.values],decoupled_setpoints.time,[decoupled_setpoints.signals.values],'--')
grid
legend('x','z','beta')
xlabel('t(s)')
ylabel('position (m), angle (rad)')
title('Decoupled Output')
subplot(2,1,2)
plot(decoupled_input.time, [decoupled_input.signals.values])
grid
legend('delta f1','delta f3','tau')
xlabel('t(s)')
ylabel('force (N), torque (Nm)')
title('Decoupled Input')
% C_theta = [0,0,1,0, 0,0,0,0];
% sim('linearized_model',te)
% sys = ss(A-B*K,V,C,D)
% impulse(sys)

