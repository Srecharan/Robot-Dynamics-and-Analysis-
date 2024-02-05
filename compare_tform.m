%% The function should return w, t, p, g_w, g_w_matlab, g_v, g_v_matlab, g and g_matlab.

function [w, t, p, g_w, g_w_matlab, g_v, g_v_matlab, g, g_matlab] = compare_tform()

    % TODO: generate a random unit vector w, rotation amount t, and the
    % twist pitch p

    w = rand(3, 1);      % Generating random unit vector w with 3x1 dimensions
    w = w/norm(w);       % Converting vector w into unit vector.
    t = 2*pi*rand();     % Generating a random rotation which lies between 0 and 2*pi
    p = 50*rand();       % Generating a random p which lies between 0 and 50

    % TODO: using w and t, compute the 4x4 rigid body transformation matrix 
    % g_w_matlab generated with axang2tform, and g_w generated with your 
    % function twist2rbvel and the Matlab function expm

    axis_angle = [w', t];                  % Represents an axis-angle representation of a rotation.
    g_w_matlab = axang2tform(axis_angle);  % g_w using twist2rbvel
    xi_w = [w*t; 0; 0; 0];                 % Angular velocity = w*t
    g_w = expm(twist2rbvel(xi_w));         % g_w using twist2rbvel

    % TODO: compute the velocity of 3-vector v = wp
    v = w*p;

    % TODO: using the pure translation with a velocity of v and the amount 
    % t, compute the 4x4 rigid body transformation matrix g_v_matlab 
    % generated with trvec2tform, and g_v generated with your function 
    % twist2rbvel and the Matlab function expm
    % Compute g_v_matlab using trvec2tform

    g_v_matlab = trvec2tform(v'*t);    % g_v_matlab using trvec2tform
    xi_v = [0; 0; 0; v*t];             % velocity = v*t
    g_v = expm(twist2rbvel(xi_v));     % g_v using twist2rbvel 

    % TODO: using w, t, and v, compute the 4x4 rigid body transformation 
    % matrix g generated with your function twist2rbvel and the Matlab 
    % function expm and g_matlab generated with the composition of 
    % axang2tform and trvec2tform
   
    g_matlab = g_w_matlab * g_v_matlab;
    xi = [w*t; v*t]; 
    g = expm(twist2rbvel(xi));

end


