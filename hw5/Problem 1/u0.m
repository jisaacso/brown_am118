function [z] = u0(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function u0 represents the    %
% I.C. of the PDE: u_t = u_xx   %
% Input:                        %
%       x = spatial location    %
% Output:                       %
%       z = solution u(x,0)     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z = sin(pi*x);