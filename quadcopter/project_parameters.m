function [mg, mq, Ig, Iq, L_center_of_mass, L_arm, gravity] = project_parameters(payload)
%% Definition of system parameters
L_arm = 0.2;
mq = 0.531;
h=0.02;
m_arm = 0.4;
L_center_of_mass=(L_arm/2*m_arm+L_arm*payload)/(m_arm+payload);
mg = m_arm+payload;
Ig = 1/12*mg*(h^2+L_arm^2)+payload*(L_arm/2)^2;
Iq = 0.00368;
% http://wiki.asctec.de/display/AR/CAD+Models
gravity = 9.81;

%% Constraints


%% Scaling factors

end