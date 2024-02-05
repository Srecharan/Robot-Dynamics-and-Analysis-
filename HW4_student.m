clear; clc; close all;

%% PROBLEM 1
% Define symbolic variables
syms q1 q2 q3 q4 q5 q6 real;
q = [q1 q2 q3 q4 q5 q6]';

% Define known frame offsets
lx2 = 320; lx5 = 887; lx6 = 200; lz2 = 680; lz3 = 975; lz4 =200;

% Compute forward  kinematics
% You can reuse your code from HW3
gs1 = [cos(q1) -sin(q1) 0 0;
       sin(q1) cos(q1) 0 0;
       0 0 1 0;
       0 0 0 1];

g12 = [cos(q2) 0 sin(q2) lx2;
       0 1 0 0;
       -sin(q2) 0 cos(q2) lz2;
       0 0 0 1];

g23 = [cos(q3) 0 sin(q3) 0;
       0 1 0 0;
       -sin(q3) 0 cos(q3) lz3;
       0 0 0 1];

g34 = [1 0 0 0;
       0 cos(q4) -sin(q4) 0;
       0 sin(q4) cos(q4) lz4;
       0 0 0 1];

g45 = [cos(q5) 0 sin(q5) lx5;
       0 1 0 0;
       -sin(q5) 0 cos(q5) 0;
       0 0 0 1];

g5t = [1 0 0 lx6;
       0 cos(q6) -sin(q6) 0;
       0 sin(q6) cos(q6) 0;
       0 0 0 1];

gst0 = [1 0 0 1407; 0 1 0 0; 0 0 1 1855; 0 0 0 1]

gst = (gs1*g12*g23*g34*g45*g5t)

%% Problem 1.1

% TODO: Calculate the inverse of Foward Kinematics transformation

gst_inv = inv(gst);

% TODO: Calculate each portion of the spatial Jacobian and get Js
J1 = sym(rbvel2twist(diff(gst,q1)*gst_inv));
J2 = sym(rbvel2twist(diff(gst,q2)*gst_inv));
J3 = sym(rbvel2twist(diff(gst,q3)*gst_inv));
J4 = sym(rbvel2twist(diff(gst,q4)*gst_inv));
J5 = sym(rbvel2twist(diff(gst,q5)*gst_inv));
J6 = sym(rbvel2twist(diff(gst,q6)*gst_inv));

Js = [J1 J2 J3 J4 J5 J6]'


%% Problem 1.2

% TODO: Calculate adjoint of product of exponential map
u1 = [0 0 1]';
u2 = [0 1 0]';
u3 = [0 1 0]';
u4 = [1 0 0]';
u5 = [0 1 0]';
u6 = [1 0 0]';

l1 = [0 0 0]';
l2 = [320 0 680]';
l3 = [320 0 1655]';
l4 = [320 0 1855]';
l5 = [1207 0 1855]';
l6 = [1407 0 1855]';

v1 = -cross(u1,l1);
v2 = -cross(u2,l2);
v3 = -cross(u3,l3);
v4 = -cross(u4,l4);
v5 = -cross(u5,l5);
v6 = -cross(u6,l6);

xi1 = [v1;u1];
xi2 = [v2;u2];
xi3 = [v3;u3];
xi4 = [v4;u4];
xi5 = [v5;u5];
xi6 = [v6;u6];

xi1_h = twist2rbvel(xi1);
xi2_h = twist2rbvel(xi2);
xi3_h = twist2rbvel(xi3);
xi4_h = twist2rbvel(xi4);
xi5_h = twist2rbvel(xi5);
xi6_h = twist2rbvel(xi6);

gst_exp = simplify(expm(xi1_h*q1)*expm(xi2_h*q2)*expm(xi3_h*q3)*expm(xi4_h*q4)*expm(xi5_h*q5)*expm(xi6_h*q6)*gst0) % Ignore

p2 = sym(expm(xi1_h*q1));
p3 = sym(expm(xi2_h*q2));
p4 = sym(expm(xi3_h*q3));
p5 = sym(expm(xi4_h*q4));
p6 = sym(expm(xi5_h*q5));

xi2_prime = sym(tform2adjoint(p2)*xi2);
xi3_prime = sym(tform2adjoint(p2*p3)*xi3);
xi4_prime = sym(tform2adjoint(p2*p3*p4)*xi4);
xi5_prime = sym(tform2adjoint(p2*p3*p4*p5)*xi5);
xi6_prime = sym(tform2adjoint(p2*p3*p4*p5*p6)*xi6);

% TODO: Calculate and compare spatial Jacobian
Js_exp = sym([xi1 xi2_prime xi3_prime xi4_prime xi5_prime xi6_prime]')

% TODO: Compare to 1.1
disp('Element-wise difference between gst and gst_exp:')
simplify(Js_exp - Js)

%% Problem 1.3

% TODO: Define body twist in initial configuration
Vb = [0 1 0 0 0 0]';

% TODO: Convert to spatial twist through adjoint, at initial configuration
Vs = tform2adjoint(gst0).*Vb

%%Problem 1.4
% TODO: Compute rank of spatial jacobian (singular if rank < 6)
Js0 = subs(Js_exp,q, [0 0 0 0 0 0]')
rank_Js0 = rank(Js0);
disp('Rank of Js0 is:')
disp(rank_Js0)

disp('In this particular configuration, the Jacobian becomes singular. This implies that rotating the 4th and 6th joints results in identical end-effector movement. Consequently, theres no specific joint velocity q dot that fulfills the equation since these joints are the only ones contributing to the y component. Furthermore, theres no combination of these or other columns that produces an exclusive y velocity without inducing motion in other directions.') 

%% PROBLEM 2
syms theta_1 theta_2 l_1 l_2 m g real

%% 2.1
% TODO: compute the forward kinematicsgstand g_st
g1 = [cos(theta_1) -sin(theta_1) 0 l_1*cos(theta_1); sin(theta_1) cos(theta_1) 0 l_1*sin(theta_1); 0 0 1 0; 0 0 0 1];
g2 = [cos(theta_2) -sin(theta_2) 0 l_2*cos(theta_2);sin(theta_2) cos(theta_2) 0 l_2*sin(theta_2); 0 0 1 0; 0 0 0 1];
g_st = simplify(g1*g2)

%% 2.2
% TODO: compute the spatial and body jacobians Js and Jb
Vbst_1 = [l_1*sin(theta_2); l_2+l_1*cos(theta_2); 0; 0; 0; 1]
Vbst_2 = [0; l_2; 0; 0; 0; 1];
J_b = [Vbst_1 Vbst_2]

gst_inv = inv(g_st);
Jx = sym(rbvel2twist(diff(g_st,theta_1)*gst_inv));
Jy = sym(rbvel2twist(diff(g_st,theta_2)*gst_inv));
J_s = [Jx Jy]'

%% 2.3
% TODO: compute the body wrench Ft as a function of the configuration
Adg_t = (tform2adjoint(g_st))';
F_s = [0; -m*g; 0; 0; 0; -m*g*(l_2*cos(theta_1 + theta_2) + l_1*cos(theta_1))];
F_t=simplify(Adg_t*F_s)


%% 2.4
% TODO: compute the joint torque tau_pm that counteracts this body wrench
tau_pm = simplify(J_b'*F_t)

%% 2.5
syms m_1 m_2 real

% TODO: compute the joint torques tau_lm that needs to be applied to counteract just the weight of the links themselves
F1 = [0;-m_1*g-m_2*g;0;0;0;-(m_1*g*(l_1/2)*cos(theta_1))-(m_2*g*((l_2/2)*cos(theta_1+theta_2)+l_1*cos(theta_1)))];
F2 = Adg_t*F1;
tau_lm = simplify(J_b'*F2)

%% 2.6
% TODO: compute the actuator effort 
energy = tau_lm'*tau_lm;
energy = simplify(energy);