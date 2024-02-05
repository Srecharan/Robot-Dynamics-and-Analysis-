%% The function should return V, g, V_s_Ad_g, V_s_tform, V_b_Ad_g and V_b_tform.

function [V, g, Vs_Adg, Vs_tform, Vb_Adg, Vb_tform] = compare_twist()
    % TODO: generate a random unit vector w, rotation amount t, and the
    % twist pitch p
    
    w = rand(3, 1);      % Generating random unit vector w with 3x1 dimensions
    w = w/norm(w);       % Converting vector w into unit vector.
    t = 2*pi*rand();     % Generating a random rotation which lies between 0 and 2*pi
    p = 50*rand();       % Generating a random p which lies between 0 and 50

    % TODO: compute the velocity of 3-vector v = wp
    v = zeros(3, 1);
    
    % TODO: using w, t, and v, construct a random twist V and compute the 
    % 4x4 rigid body transformation matrix g
    
    axis_angle = [w', t];
    g_w = axang2tform(axis_angle);
    v = w * p;
    g_v = trvec2tform(v'*t);
    g = g_w * g_v;

    V = rand(6, 1);     % Generating random unit vector w with 6x1 dimensions

     % TODO: first, treating V as a body velocity, compute the conversion to 
    % spatial velocity V_s_Ad_g using tform2adjoint and the conversion in 
    % homogeneous coordinates V_s_tform using twist2rbvel and rbvel2twist
    
    Adg = tform2adjoint(g);    %  Body velocity into spatial velocity using Adjoint transformation
    Vs_Adg = Adg * V;

    Vs_tform = rbvel2twist(twist2rbvel(Vs_Adg));

    Adg_inv = inv(Adg);        %  Spatial velocity to body velocity using Adjoint transformation
    Vb_Adg = Adg_inv * V;

    Vb_tform = rbvel2twist(twist2rbvel(Vb_Adg));

end


