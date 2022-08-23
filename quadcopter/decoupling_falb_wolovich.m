clear
symbolic_rotor()
B=B_dbth;
% B(5,3)=0;
n=8;
p=3;
delta_1 = zeros(p,1);
for i=1:p
    for j=0:100
        delta_1(i)=j;
        ci = C(i,:);
        judge = (ci*A^delta_1(i)*B==zeros(1,3));
        if any(judge<0.5)
            break;
        end
    end
end
delta = delta_1+1;
D_star0 = zeros(p,p);
for i=1:p
    ci = C(i,:);
    D_star0(i,:)=ci*A^delta_1(i)*B;
end
%although y1=x has the order 2 because of B53, its major order is
%4 due to the indirect control of M-theta_dot-theta(x_dotdot)-x_dot-x
delta_1(1)=3;
delta(1)=4;
D_star = zeros(p,p);
for i=1:p
    ci = C(i,:);
    D_star(i,:)=ci*A^delta_1(i)*B;
end
K=diag([512,9,16]);
F=[C(1,:)*A^delta(1)+512*C(1,:)*A^0+704*C(1,:)*A+216*C(1,:)*A^2+25*C(1,:)*A^3;
    C(2,:)*A^delta(2)+9*C(2,:)*A^0+6*C(2,:)*A;
    C(3,:)*A^delta(3)+16*C(3,:)*A^0+8*C(3,:)*A];
R=D_star\F;
M=D_star\K;
x0=zeros(8,1);
te=14;
sim('linearized_decoupling',te)