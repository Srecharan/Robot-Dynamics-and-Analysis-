function dA = computedA(in1,in2)
%COMPUTEDA
%    DA = COMPUTEDA(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    04-Nov-2019 15:37:45

dx = in2(1,:);
dy = in2(2,:);
dA = reshape([0.0,0.0,dx.*2.0,0.0,0.0,dy.*2.0],[3,2]);