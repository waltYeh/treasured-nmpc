function plot_rings(radius,yaw,x_tr,y_tr,z_tr)
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
plot3(x_stem,y_stem,z_stem,'ko')
axis equal
end