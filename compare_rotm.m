%% The function should return w, t, rotm and rotm_matlab, and rotm and rotm_matlab should be the same.

function [w, t, rotm, rotm_matlab] = compare_rotm()

    % TODO: generate a random unit vector w and rotation amount t

    w = rand(3, 1);   % Generating random unit vector w with 3x1 dimensions
    w = w/norm(w);    % Converting vector w into unit vector
    t = 2*pi*rand();  % Generating a random rotation which lies between 0 and 2*pi

    % TODO: Compute the 3x3 rotation matrix rotm_matlab generated with 
    % Matlab built-in function axang2rotm 

    axis_angle = [w', t];                  % Represents an axis-angle representation of a rotation
    rotm_matlab = axang2rotm(axis_angle);  % Axis angle is made into 3x3 matrix and stored in rotation matrix

    % TODO: Compute the 3x3 rotation matrix rotm generated with generated 
    % with your function angvel2skew and the Matlab function expm

    skew_mat = angvel2skew(w * t);     % Angular velocity = w*t
    rotm = expm(skew_mat);             

end
