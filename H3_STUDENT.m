clear all; clc; close all;
%% Problem 2.1
% TODO: Writing out by hand the rigid body transformation from the stationary 
% frame to the tool frame in the initial configuration

gst0 = [1, 0, 0, 1407;  % Rigid body transformation from the figure
       0, 1, 0, 0;
       0, 0, 1, 1855;
       0, 0, 0, 1];
%% Problem 2.2

% Define symbolic variables
syms q1 q2 q3 q4 q5 q6 real;
q = [q1 q2 q3 q4 q5 q6]';

% Define known frame offsets
lx2 = 320; lx5 = 887; lx6 = 200; lz2 = 680; lz3 = 975; lz4 =200;

% TODO: Define rigid body transformations between successive links
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

% TODO: Compute forward kinematics
gst = simplify(gs1*g12*g23*g34*g45*g5t); % Calculating the general forward kinematics matrix 

% TODO: Compute gst(0) with assigning all symbolic variables as 0
gst0_sym = double(subs(gst, q, zeros(6, 1))); % gst0 is the matrix as found in part 2.1

% TODO: Compare gst0 to gst0_sym
disp('Element-wise difference between gst0 and gst0_sym: ')
disp(gst0 - gst0_sym)
%% Problem 2.3

% TODO: Define joint twists in initial configuration
xi1 = [0;
       0;
       0;
       0;
       0;
       1];

xi2 = [-lz2;
       0;
       lx2;
       0;
       1;
       0];

xi3 = [-(lz2+lz3);
       0;
       lx2;
       0;
       1;
       0];

xi4 = [0;
       lz2+lz3+lz4;
       0;
       1;
       0;
       0];

xi5 = [-(lz2+lz3+lz4);
       0;
       lx2+lx5;
       0;
       1;
       0];

xi6 = [0;
       lz2+lz3+lz4;
       0;
       1;
       0;
       0];

% TODO: Compute product of exponentials
gst_exp = expm(twist2rbvel(xi1)*q1)*expm(twist2rbvel(xi2)*q2)*expm(twist2rbvel(xi3)*q3)*expm(twist2rbvel(xi4)*q4)*expm(twist2rbvel(xi5)*q5)*expm(twist2rbvel(xi6)*q6)*gst0;
gst_exp = simplify(gst_exp);

% TODO: Compare to 2.2
disp('Element-wise difference between gst and gst_exp:')
disp(simplify(gst_exp - gst))
%% Problem 2.4

% TODO: Define goal position as initial configuration +100mm in +y direction
gDes = gst0;
gDes(2, 4) = gDes(2, 4) + 100;

% TODO: Define an optimization cost function
F = sum(sum((gst-gDes).^2));

% Define an optimization problem
fun = matlabFunction(F, 'var', {q});

% Call fminunc to solve Inverse Kinematics
options = optimoptions(@fminunc, 'Display', 'iter');
q_sol = fminunc(fun, zeros(6, 1), options);

% TODO: Compute the Forward Kinematics using the IK solution
gAchieved = double(subs(gst, q, q_sol));

% TODO: Compare gAchieved with gDes
disp('Element-wise difference between gst and gst_exp:')
disp(gAchieved - gDes)