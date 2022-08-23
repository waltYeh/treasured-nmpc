y=zeros(2*N-2,1);
Phi=zeros(2*N-2,3);
for k = 4:2*N+1
    y(k-3)=w(k)+w(k-2);
    Phi(k-3,1) = M(k)-M(k-2);
    Phi(k-3,2) = acc(k)-acc(k-2);
    Phi(k-3,3) = -w(k-1);
end
theta0=Phi\y;
alpha_ref = 2*T/(4*I+T^2*mg*gravity*L);
theta_ref = [alpha_ref;mg*L*alpha_ref;(-8*I+2*T^2*mg*gravity*L)/(4*I+T^2*mg*gravity*L)];


yalg = [0;0];
Malg=zeros(2,2);
Malg(1,1) = 0;
int_a = zeros(2*N,1);
int_omg=zeros(2*N,1);
intint_omg=zeros(2*N,1);
intint_a=0;
intint_M=0;
intintint_omg = 0;
int_M=zeros(2*N,1);
for i=1:2*N
    int_a(i+1) = int_a(i)+T*acc(i);
    int_omg(i+1) = int_omg(i)+T*w(i);
    int_M(i+1) = int_M(i) +T*M(i);
end
for i=1:2*N
    intint_omg(i+1) = intint_omg(i)+T*int_omg(i);
    intint_a = intint_a+T*int_a(i);
    intint_M = intint_M+T*int_M(i);
end
for i=1:2*N
    intintint_omg = intintint_omg+T*intint_omg(i);
end
Malg(1,2) = -(int_a(end)+gravity*intint_omg(end));
yalg(1) = int_M(end);
Malg(2,1) = int_omg(end);
Malg(2,2) = -(intint_a(end)+gravity*intintint_omg(end));
yalg(2) = intint_M(end);
theta_alg = Malg\yalg;
theta_alg_ref = [I;mg*L];