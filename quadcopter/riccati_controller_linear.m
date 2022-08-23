clear
symbolic_rotor()
x0=zeros(8,1);
Q=diag([5,25,0.1,20, 1,3,1,1]);
% more for x,z: faster towards desired pos, but with overshoot
% more for x_dot, z_dot; slower with no overshoot
% more for beta: beta less affected by positioning
% more for beta_dot: beta tremble less
% more for theta_dot: less difference of motors, x slower
R=diag([100,100,10]);
% more f1, f3: slow down
% more tau: worse beta for lowering cost of tau, less peak of tau
B_here = B_dbth;
%% calculation of the solution of the riccati equation and the controller

[P,L,K] = care(A,B_here,Q,R,zeros(size(B)),eye(size(A)));
V=inv(C*((B_here*K-A)\B_here));
te=10;
sim('linearized_coupled_model',te)
figure(1)
subplot(2,1,1)
plot(coupled_output.time, [coupled_output.signals.values],coupled_setpoints.time,[coupled_setpoints.signals.values],'--')
grid
legend('x','z','beta')
xlabel('t(s)')
ylabel('position (m), angle (rad)')
title('Coupled Output')
subplot(2,1,2)
plot(coupled_input.time, [coupled_input.signals.values])
grid
legend('delta f1','delta f3','tau')
xlabel('t(s)')
ylabel('force (N), torque (Nm)')
title('Coupled Input')
% sys = ss(A-B*K,V,C,D)
% impulse(sys)
