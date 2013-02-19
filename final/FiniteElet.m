function U = FiniteElet(x)

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FinitElet uses a Finite Elmement method      %
% to solve any ODE of the form u(t)''=f(t)     %
% whose Right Hand Side, f(t) is stored        %
% within the .m file "f" and with zero B.C.'s  %
% using uniform grid sizing                    %
% INPUT                                        %
%    x  = partition of space                   %
%          x(1)=x0,x(2)=x1,...,x(N+1)=xN       %
%                                              %
%                                              %
% OUTPUT:                                      %
%    err = error of approximation at time T    %
% @author: Joseph Isaacson
% -Brown Univerisy Undergraduate
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %first, error check:

 if(x(1)~=0 || x(length(x))~=1)
     ERROR = 'Error! X(0) must = 0 and X(N) must =1!!!'
     return
 end
 
 %initialize variables:
 N = length(x)-1;      %N = number of time steps (10)
 U = zeros(N-2,1);    %U(x) is the approx solution at spatial step x
 
 %construct A:
 for(i=1:N-1)
     c1(i) = (x(i+2)-x(i))/((x(i)-x(i+1))*(x(i+1)-x(i+2)));
 end
 for(i=1:N-2)
    c2(i) = 1/(x(i+1)-x(i+2));
 end
 A = diag(c1) + diag(c2,1) + diag(c2,-1);
 A = sparse(A);
 
 %construct F
 F = zeros(N-1,1);
 for(i = 1:N-1)
%      midpoint method to approx (f,phi(i))
     c1 = x(i+1)-x(i);
     c2 = (x(i+1)+x(i))/2;
     c3 = x(i+2)-x(i+1);
     c4 = (x(i+2)+x(i+1))/2;
     F(i) = c1*f(c2)*phi(i,c2,x)+c3*f(c4)*phi(i,c4,x);
 end
 
 U = (A\F);
 U = [0,U',0];
 
 hold on;
 plot(x,U,'k');
 plot([0:1/10000:1],sin(pi*[0:1/10000:1]))
 hold off;
 
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function u0 represents the    %
% I.C. of the PDE: u_t = u_xx   %
% Input:                        %
%       x = spatial location    %
% Output:                       %
%       z = solution u(x,0)     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 function [z] = f(x)
 z = pi^2*sin(pi*x);



%part = partition of space (= x in FiniteElet.m)
%x = current spatial location
%j = current phi
function z = phi(j,x,part)

if(part(j)<=x && part(j+1)>=x)
    z = (x-part(j))./(part(j+1)-part(j));
elseif(x>part(j+1) && x<part(j+2))
    z = (part(j+2)-x)./(part(j+2)-part(j+1));
else
    z = 0;
end


function z = integrand(j,x,part)
z=(f(x).*phi(j,x,part));