function [z] = exact(x,t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function exact represents the %
% exact solution of u_t = u_xx  %
% with zero B.C.'s and with an  %
% I.C. stored in u0(x)          %
% Input:                        %
%       x = spatial location    %
%       t = location in time    %
% Output:                       %
%       z = solution u(x,t)     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z = 0;
for(j = 1:10)   %use 10 terms in fourier series
    c = 2*quad(@(x)integrand(x,j),0,1);
    z = z + c*sin(j*pi*x)*exp(-j^2*pi^2*t);
end