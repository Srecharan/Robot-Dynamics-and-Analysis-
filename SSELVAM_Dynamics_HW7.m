%% Problem 1 Initialization

clear;
clc;

syms q1 q2 dq1 dq2 ddq1 ddq2 real
syms m g l real
syms F tau real
q_g = [q1;q2];
dq_g = [dq1;dq2];
ddq_g = [ddq1;ddq2];
I_o = 1/12*m*l^2;

%% Problem 1.1

v1 = dq1*l/2;
v2 = sqrt((dq1*(l/2 + q2))^2 + dq2^2);

% TODO: Kinematic energy
T_g = simplify(1/2*m*v1^2 + 1/2*m*v2^2 + 1/2*I_o*dq1^2 + 1/2*I_o*dq1^2);

% TODO: Potential energy
V_g = simplify(m*g*l/2*sin(q1) + m*g*(l/2 + q2)*sin(q1));

% TODO: Lagrangian in generalized coordinates
L_g = T_g-V_g;

%% Problem 1.2

% TODO: Applied forces to each generalized coordinate
Y_g = sym(zeros(2, 1));
Y_g = [tau;F];

% TODO: Equations of motion by using the Lagrange equations to differentiate the Lagrangian in generalized coordinates
EOM_g1 = sym(zeros(2, 1));
EOM_g1 = simplify(jacobian(jacobian(L_g,dq_g)',[q_g;dq_g])*[dq_g;ddq_g]-jacobian(L_g,q_g)'-Y_g);

%% Problem 1.3

M1 = [m 0 0;
      0 m 0;
      0 0 I_o];
M2 = M1;

% For computing the kinetic energy, we need to do
Jbsl1 = [0 0;
          l/2 0
          1 0];
Jbsl2 = [0 1
          l/2 + q2 0
          1 0];

% TODO: Inertia tensor
M_g = Jbsl1'*M1*Jbsl1 + Jbsl2'*M2*Jbsl2;

% TODO: Coriolis matrix
C_g  = sym(zeros(length(q_g),length(q_g)));
for i = 1:length(q_g)
    for j = 1:length(q_g)
        for k = 1:length(q_g)
            C_g(i,j) = C_g(i,j) + 1/2*(diff(M_g(i,j),q_g(k)) + diff(M_g(i,k),q_g(j)) - diff(M_g(j,k),q_g(i)))*dq_g(k);
        end
    end
end
C_g = simplify(C_g);

% TODO: Nonlinear terms
N_g = [diff(V_g, q1);diff(V_g, q2)];

% TODO: Equations of motion computed directly
EOM_g2 = M_g*ddq_g + C_g*dq_g + N_g - Y_g;

%% Problem 2 Initialization

% Declare symbolic variables for maximal coordinates and collect into vectors
syms x1 y1 phi_1 x2 y2 phi_2 dx1 dy1 dphi_1 dx2 dy2 dphi_2 ddx1 ddy1 ddphi_1 ddx2 ddy2 ddphi_2 real
syms lambda_1 lambda_2 lambda_3 lambda_4 real
q_m = [x1; y1; phi_1; x2; y2; phi_2];
dq_m = [dx1; dy1; dphi_1; dx2; dy2; dphi_2];
ddq_m = [ddx1; ddy1; ddphi_1; ddx2; ddy2; ddphi_2];
lambda = [lambda_1; lambda_2; lambda_3; lambda_4];
I_o = 1/12*m*l^2;

%% Problem 2.1

% TODO: Constraint function
a_m = sym(zeros(4, 1));
a_m = [x1 - l/2*cos(phi_1)
    y1 - l/2*sin(phi_1)
    y2*cos(phi_2) - x2*sin(phi_2)
    phi_2 - phi_1];

% TODO: A as the differential of a
A_m = sym(zeros(4, 6));
A_m = jacobian(a_m, q_m);

%% Problem 2.2

% TODO: Kinematic energy
T_m = sym(0);
T_m = 1/2*m*(dx1^2 + dy1^2 + dx2^2 + dy2^2) + 1/2*I_o*(dphi_1^2 + dphi_2^2);

% TODO: Potential energy
V_m = sym(0);
V_m = m*g*(y1 + y2);

% TODO: Lagrangian in maximal coordinates
L_m = sym(0);
L_m = T_m-V_m;

%% Problem 2.3

% TODO: Applied forces in maximal coordinates
Y_m = sym(zeros(6, 1));
Y_m = [-F*cos(phi_1); -F*sin(phi_1); tau; F*cos(phi_1); F*sin(phi_1); 0];

%% Problem 2.4

% TODO: Equations of motion by using the Lagrange equations to differentiate the Lagrangian in maximal coordinates
EOM_m1 = sym(zeros(6, 1));
EOM_m1 = simplify(jacobian(jacobian(L_m,dq_m)',[q_m;dq_m])*[dq_m;ddq_m]-jacobian(L_m,q_m)'+A_m'*lambda-Y_m);

%% Problem 2.5

M1 = [m 0 0;0 m 0;0 0 I_o];
M2 = M1;

% Compute body Jacobians for each link. See PDF solution for reasoning.
Jbsl1 = [cos(phi_1) sin(phi_1) 0 0 0 0;
        -sin(phi_1) cos(phi_1) 0 0 0 0
         0 0 1 0 0 0];
Jbsl2 = [0 0 0 cos(phi_2) sin(phi_2) 0;
         0 0 0 -sin(phi_2) cos(phi_2) 0
         0 0 0 0 0 1];

% TODO: Inertia tensor
M_m = sym(zeros(6, 6));
M_m = simplify(Jbsl1'*M1*Jbsl1 + Jbsl2'*M2*Jbsl2);

% TODO: Coriolis matrix
C_m  = sym(zeros(6, 6));
 for i = 1:length(q_m)
     for j = 1:length(q_m)
         for k = 1:length(q_m)
             C_m(i,j) = C_m(i,j) + 1/2*(diff(M_m(i,j),q_m(k)) + diff(M_m(i,k),q_m(j)) - diff(M_m(j,k),q_m(i)))*dq_m(k);
         end
     end
 end

% TODO: Nonlinear terms
N_m = sym(zeros(6, 1));
N_m = jacobian(V_m,q_m)';

% TODO: Equations of motion computed directly
EOM_m2 = sym(zeros(6, 1));
EOM_m2 = simplify(M_m*ddq_m + C_m*dq_m + N_m  + A_m'*lambda - Y_m);

%% Problem 2.6

% TODO: Constraint forces
lambdaVec = sym(zeros(4, 1));
dA = sym(zeros(4,6));
for i = 1:length(q_m)
    dA = dA + diff(A_m, q_m(i))*dq_m(i);
end

Minv = inv(M_m);
lambdaVec = simplify((A_m*Minv*A_m')\(A_m*Minv*(Y_m - C_m*dq_m - N_m) + dA*dq_m));

%% Problem 2.7 (Optional)

% TODO: q_m = h(q_g)
h = sym(zeros(6, 1));

% TODO: dq_m = H*dq_g
H = sym(zeros(6, 2));

% TODO: Reconstruct
Mhat = sym(zeros(2, 2));
Chat = sym(zeros(2, 2));
Nhat = sym(zeros(2, 1));
Yhat = sym(zeros(2, 1));

% TODO: Equations of motion
EOM_g3 = sym(zeros(2, 1));
