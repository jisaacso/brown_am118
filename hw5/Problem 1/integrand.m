function [z] = integrand(x,j)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function integrand represents  %
% the integrand when finding the %
% fourier coefficients of the    %
% exact solution of u_t = u_xx   %
% with zero B.C.'s and with an   %
% I.C. stored in u0(x)           %
% Input:                         %
%       x = spatial location     %
%       j = current term in      %
%           fourier series       %
% Output:                        %
%       z = integrand            %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z = u0(x).*sin(j*pi*x);