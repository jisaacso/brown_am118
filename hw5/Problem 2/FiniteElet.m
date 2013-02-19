function err = FiniteElet(t0,T,h)

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FinitElet uses a Finite Elmement method      %
% to solve any ODE of the form u(t)''=f(t)     %
% whose Right Hand Side, f(t) is stored        %
% within the .m file "f" and with zero B.C.'s  %
% using uniform grid sizing                    %
% INPUTS:                                      %
%    t0 = initial time                         %
%    T  = final time                           %
%    h  = step size                            %
%                                              %
%                                              %
% OUTPUT:                                      %
%    err = error of approximation at time T    %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
 %initialize variables:
 N = (T-t0)/h;      %N = number of time steps (10)
 x = linspace(t0,T,N+1); % x = partition of space (0,.1,.2,.3,.4,...,1) (11)
 U = zeros(N-2,1);    %U(x) is the approx solution at spatial step x
 
 %construct A:
 %Note: since we are taking uniform interval sizes, the areas of all
 %of the basis functions will be equal, hence all diagonal values of A
 %will be equal (as will the off diagonals)

 c1 = (x(3)-x(1))/(x(1)-x(2))^2;
 c2 = 1/(x(2)-x(3));
 A = diag(c1*ones(N-1,1)) + diag(c2*ones(N-2,1),1) +...
        diag(c2*ones(N-2,1),-1);
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

 err = 0;
 %calculate error:
 for(i=1:N)
    err = err + abs(sin(pi/2*(x(i)+x(i+1)))-( (U(i)+U(i+1)) /2 ) )^2;
 end
 err = sqrt(h*err);
 
 hold on;
 plot(x,U,'k');
 plot([0:1/10000:1],sin(pi*[0:1/10000:1]))
 hold off;