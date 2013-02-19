function TolErr = FD_Explicit_PDE(t0,T,h,r)

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
 U = zeros(m,1);    %U(x) is the approx solution at spatial step x
 
 %initialize U0
 for(i = 1:m) 
     U(i) = u0((i-1)*h);
 end
 
 
 %construct A
 A = sparse(m,m);
 A(1,1) = 1-2*r;
 A(1,2) = r;
 A(m,m) = 1-2*r;
 A(m,m-1) = r;
 for(i = 2:m-1)
     A(i,i) = 1-2*r;
     A(i,i+1) = r;
     A(i,i-1) = r;
 end
 
 %find U
 for(i = 1:N)
     U = A*U;
     U(1) = 0;
     U(m) = 0;
 end
 
 %calculate error at time = T
 err = 0;
 for(j = 1:m)
    u = exact(j*h,T);           %find exact soln: u(x,T)
    err = err + abs(u-U(j))^2;  %sum up error using norm from CH2.
 end

 TolErr = sqrt(h)*sqrt(err);    %output total error
 
     
     
     
     
     
 
 