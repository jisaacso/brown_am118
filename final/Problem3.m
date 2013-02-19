function TolErr = Problem3(t0,T,h,r)

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FD_Implicit uses an Explicit finite difference %
% method (Euler in time) to solve the PDE:       %
% u_t = u_xx with u(0,t) = u(1,t) = 0            %
% and where u(x,0) = f(x) is stored within the   %
% .m file "u0"                                   %
% INPUTS:                                        %
%    t0 = initial time                           %
%    T  = final time                             %
%    h  = spatial step size                      %
%    r  = ratio of time to spatial steps         %
%         (from CFL condition) r = k/h^2         %
%                                                %
%    sample command line call:                   %
%    FD_Explicit_PDE(0,1,1/64,.4)                %
%                                                %
% OUTPUT:                                        %
%    TolErr = Error of FD method at time =  T    %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 %initialize variables:
 k = r*h^2;         %k = time step size
 m = 1/h;           %m = number of spatial steps
 N = (T-t0)/k;      %N = number of time steps
 U = zeros(m+1,N+1);  %U(x,t) is the approx solution at spatial step x
                    %and time t   
 
 %initialize U0
 for(i = 2:m) 
     U(i,1) = u0((i-1)*h);
 end
 
 
 %construct A
 A = sparse(m+1,m+1);
 for(i = 2:m)
     A(i,i) = 1-2*r;
     A(i,i+1) = r;
     A(i,i-1) = r;
 end
 
 %find U at t = 1
 U(:,2) = A*U(:,1);
 
 %reconstruct A
 A = sparse(m+1,m+1);
 for(i = 2:m)
     A(i,i) = -4*r;
     A(i,i+1) = 2*r;
     A(i,i-1) = 2*r;
 end
 
 %find U for all t=2 to t=N
 hold on;
 for(i=1:N-1)
     plot(U(:,i));
     drawnow();
     U(:,i+2)=U(:,i)+A*U(:,i+1);
     pause();
 end
 hold off;
 
 
 
 %calculate error at time = T
 err = 0;
 for(j = 1:m)
    u = exact(j*h,T);           %find exact soln: u(x,T)
    err = err + abs(u-U(j,N+1))^2;  %sum up error using norm from CH2.
 end

 TolErr = sqrt(h)*sqrt(err);    %output total error
 
 
 
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
     
     
 
 