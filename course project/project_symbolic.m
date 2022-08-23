clear
syms t
syms q1(t) q2(t)
syms a1 a2 a3
syms bb1 bb2 bb3
syms c1 c2 s1 s2
syms l1 l2
syms m1 m2 m3 g
syms J1 J2
syms u1 u2
syms q1_dot q2_dot q1_ddot q2_ddot

b1=cos(q1)*a1+sin(q1)*a2;
b2=-sin(q1)*a1+cos(q1)*a2;
b3=a3;
% c1=cos(q2)*b1+sin(q2)*b2;
% c2=-sin(q2)*b1+cos(q2)*b2;
% c3=b3;
% q1_dot=diff(q1,t);
% q2_dot=diff(q2,t);
w_A_B=q1_dot*a3;
w_A_C=(q1_dot+q2_dot)*a3;
r_A_P=l1/2*b1;
r_A_Q=l1*b1+l2/2*(cos(q2)*b1+sin(q2)*b2);
r_A_R=l1*b1+l2*(cos(q2)*b1+sin(q2)*b2);
v_A_P=diff(r_A_P,t);
v_A_Q=diff(r_A_Q,t);
v_A_R=diff(r_A_R,t);
a_A_P=diff(v_A_P,t);
a_A_Q=diff(v_A_Q,t);
a_A_R=diff(v_A_R,t);




dw_A_P_dq1_dot=diff(subs(w_A_B,diff(q1(t), t),q1_dot),q1_dot);
dw_A_P_dq1_dot=subs(dw_A_P_dq1_dot,a1,cos(q1)*bb1-sin(q1)*bb2);
dw_A_P_dq1_dot=subs(dw_A_P_dq1_dot,a2,sin(q1)*bb1+cos(q1)*bb2);
dw_A_P_dq1_dot=subs(dw_A_P_dq1_dot,a3,bb3);
dw_A_P_dq1_dot=simplify(dw_A_P_dq1_dot);
dw_A_P_dq1_dot=collect(dw_A_P_dq1_dot,{'bb1','bb2','bb3'});
dw_A_P_dq2_dot=0;
dw_A_Q_dq1_dot=bb3;
dw_A_Q_dq2_dot=bb3;



% dv_A_P_dq1dot=-l1/2*s1*a1+l1/2*c1*a2;
dv_A_P_dq1_dot=diff(subs(v_A_P,diff(q1(t), t),q1_dot),q1_dot);
dv_A_P_dq1_dot=subs(dv_A_P_dq1_dot,a1,cos(q1)*bb1-sin(q1)*bb2);
dv_A_P_dq1_dot=subs(dv_A_P_dq1_dot,a2,sin(q1)*bb1+cos(q1)*bb2);
dv_A_P_dq1_dot=subs(dv_A_P_dq1_dot,a3,bb3);
dv_A_P_dq1_dot=simplify(dv_A_P_dq1_dot);
dv_A_P_dq1_dot=collect(dv_A_P_dq1_dot,{'bb1','bb2','bb3'});

dv_A_P_dq2_dot=0;
% dv_A_Q_dq1dot=-l2/2*s2*b1+(l1+l2/2*c2)*b2;
dv_A_Q_dq1_dot=diff(subs(v_A_Q,diff(q1(t), t),q1_dot),q1_dot);
dv_A_Q_dq1_dot=subs(dv_A_Q_dq1_dot,a1,cos(q1)*bb1-sin(q1)*bb2);
dv_A_Q_dq1_dot=subs(dv_A_Q_dq1_dot,a2,sin(q1)*bb1+cos(q1)*bb2);
dv_A_Q_dq1_dot=subs(dv_A_Q_dq1_dot,a3,bb3);
dv_A_Q_dq1_dot=simplify(dv_A_Q_dq1_dot);
dv_A_Q_dq1_dot=collect(dv_A_Q_dq1_dot,{'bb1','bb2','bb3'});

% dv_A_Q_dq2dot=-l2/2*s2*b1+l2/2*c2*b2;
dv_A_Q_dq2_dot=diff(subs(v_A_Q,diff(q2(t), t),q2_dot),q2_dot);
dv_A_Q_dq2_dot=subs(dv_A_Q_dq2_dot,a1,cos(q1)*bb1-sin(q1)*bb2);
dv_A_Q_dq2_dot=subs(dv_A_Q_dq2_dot,a2,sin(q1)*bb1+cos(q1)*bb2);
dv_A_Q_dq2_dot=subs(dv_A_Q_dq2_dot,a3,bb3);
dv_A_Q_dq2_dot=simplify(dv_A_Q_dq2_dot);
dv_A_Q_dq2_dot=collect(dv_A_Q_dq2_dot,{'bb1','bb2','bb3'});

% dv_A_R_dq1dot=-l2*s2*b1+(l1+l2*c2)*b2;
dv_A_R_dq1_dot=diff(subs(v_A_R,diff(q1(t), t),q1_dot),q1_dot);
dv_A_R_dq1_dot=subs(dv_A_R_dq1_dot,a1,cos(q1)*bb1-sin(q1)*bb2);
dv_A_R_dq1_dot=subs(dv_A_R_dq1_dot,a2,sin(q1)*bb1+cos(q1)*bb2);
dv_A_R_dq1_dot=subs(dv_A_R_dq1_dot,a3,bb3);
dv_A_R_dq1_dot=simplify(dv_A_R_dq1_dot);
dv_A_R_dq1_dot=collect(dv_A_R_dq1_dot,{'bb1','bb2','bb3'});

% dv_A_R_dq2dot=-l2*s2*b1+l2*c2*b2;
dv_A_R_dq2_dot=diff(subs(v_A_R,diff(q2(t), t),q2_dot),q2_dot);
dv_A_R_dq2_dot=subs(dv_A_R_dq2_dot,a1,cos(q1)*bb1-sin(q1)*bb2);
dv_A_R_dq2_dot=subs(dv_A_R_dq2_dot,a2,sin(q1)*bb1+cos(q1)*bb2);
dv_A_R_dq2_dot=subs(dv_A_R_dq2_dot,a3,bb3);
dv_A_R_dq2_dot=simplify(dv_A_R_dq2_dot);
dv_A_R_dq2_dot=collect(dv_A_R_dq2_dot,{'bb1','bb2','bb3'});

Force_P=m1*a_A_P+m1*g*a2;
Force_P=subs(Force_P,a1,cos(q1)*bb1-sin(q1)*bb2);
Force_P=subs(Force_P,a2,sin(q1)*bb1+cos(q1)*bb2);
Force_P=subs(Force_P,a3,bb3);
Force_P=simplify(Force_P);
Force_P=collect(Force_P,{'bb1','bb2','bb3'});

Force_Q=m2*a_A_Q+m2*g*a2;
Force_Q=subs(Force_Q,a1,cos(q1)*bb1-sin(q1)*bb2);
Force_Q=subs(Force_Q,a2,sin(q1)*bb1+cos(q1)*bb2);
Force_Q=subs(Force_Q,a3,bb3);
Force_Q=simplify(Force_Q);
Force_Q=collect(Force_Q,{'bb1','bb2','bb3'});

Force_R=m3*a_A_R+m3*g*a2;
Force_R=subs(Force_R,a1,cos(q1)*bb1-sin(q1)*bb2);
Force_R=subs(Force_R,a2,sin(q1)*bb1+cos(q1)*bb2);
Force_R=subs(Force_R,a3,bb3);
Force_R=simplify(Force_R);
Force_R=collect(Force_R,{'bb1','bb2','bb3'});

H_P=diff(J1*q1,t)*bb3;
H_Q=diff(J2*(q1+q2),t)*bb3;
shall_zero2=diff(Force_P,bb1)*diff(dv_A_P_dq2_dot,bb1)+diff(Force_P,bb2)*diff(dv_A_P_dq2_dot,bb2)+diff((diff(H_P)-(u1-u2)*bb3),bb3)*diff(dw_A_P_dq2_dot,bb3)...%arm1
    +diff(Force_Q,bb1)*diff(dv_A_Q_dq2_dot,bb1)+diff(Force_Q,bb2)*diff(dv_A_Q_dq2_dot,bb2)+diff((diff(H_Q)-u2*bb3),bb3)*diff(dw_A_Q_dq2_dot,bb3)...%arm2
    +diff(Force_R,bb1)*diff(dv_A_R_dq2_dot,bb1)+diff(Force_R,bb2)*diff(dv_A_R_dq2_dot,bb2);%end
shall_zero2=simplify(shall_zero2);
shall_zero2=subs(shall_zero2,diff(q1(t), t),q1_dot);
shall_zero2=subs(shall_zero2,diff(q2(t), t),q2_dot);
shall_zero2=subs(shall_zero2,diff(q1(t), t, t),q1_ddot);
shall_zero2=subs(shall_zero2,diff(q2(t), t, t),q2_ddot);
shall_zero2=subs(shall_zero2,cos(q1(t)),c1);
shall_zero2=subs(shall_zero2,sin(q1(t)),s1);
shall_zero2=subs(shall_zero2,cos(q2(t)),c2);
shall_zero2=subs(shall_zero2,sin(q2(t)),s2);
shall_zero2=collect(shall_zero2,{'q1_ddot','q2_ddot','q1_dot','q2_dot'});

shall_zero1=diff(Force_P,bb1)*diff(dv_A_P_dq1_dot,bb1)+diff(Force_P,bb2)*diff(dv_A_P_dq1_dot,bb2)+diff((diff(H_P)-(u1-u2)*bb3),bb3)*diff(dw_A_P_dq1_dot,bb3)...%arm1
    +diff(Force_Q,bb1)*diff(dv_A_Q_dq1_dot,bb1)+diff(Force_Q,bb2)*diff(dv_A_Q_dq1_dot,bb2)+diff((diff(H_Q)-u2*bb3),bb3)*diff(dw_A_Q_dq1_dot,bb3)...%arm2
    +diff(Force_R,bb1)*diff(dv_A_R_dq1_dot,bb1)+diff(Force_R,bb2)*diff(dv_A_R_dq1_dot,bb2);%end
shall_zero1=simplify(shall_zero1);
shall_zero1=subs(shall_zero1,diff(q1(t), t),q1_dot);
shall_zero1=subs(shall_zero1,diff(q2(t), t),q2_dot);
shall_zero1=subs(shall_zero1,diff(q1(t), t, t),q1_ddot);
shall_zero1=subs(shall_zero1,diff(q2(t), t, t),q2_ddot);
shall_zero1=subs(shall_zero1,cos(q1(t)),c1);
shall_zero1=subs(shall_zero1,sin(q1(t)),s1);
shall_zero1=subs(shall_zero1,cos(q2(t)),c2);
shall_zero1=subs(shall_zero1,sin(q2(t)),s2);
shall_zero1=collect(shall_zero1,{'q1_ddot','q2_ddot','q1_dot','q2_dot'});
