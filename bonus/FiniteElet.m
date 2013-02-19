function err = FiniteElet(X)

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
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %first, error check:

 if(X(1)~=0 || X(length(X))~=1)
     error('X(0) must = 0 and X(N) must =1!!!');
 end
 
 
 %initialize variables:
 N = 2*length(X)-3;      %N = number of basis functions
 U = zeros(N-2,1);    %U(x) is the approx solution at spatial step x
 Y = (X(1:length(X)-1)+X(2:length(X)))/2; %Y is midpoint of intervals
 tol = 100*eps; %integral tolerances
 
 
 %construct A:
 for(i=1:N) 
     %calculate (phi(i),phi(i))
     odd = (i+1)/2-1; 
     even = i/2-1;
     if(mod(i,2)~=0) %if odd phi
         a = X(1+odd);
         b = X(2+odd);
     else
         a = X(1+even);
         b = X(3+even);
     end
     c1(i) = quadv(@(x)(phideriv(i,x,X,Y))^2,a,b,tol);
 end
 for(i=1:N-1) 
     %calculate (phi(i),phi(i+1))
     odd = (i+1)/2-1; 
     even = i/2-1;
     if(mod(i,2)~=0) %if odd phi
         a = X(1+odd);
         b = X(2+odd);
     else
         a = X(2+even);
         b = X(3+even);
     end
     c2(i) = quadv(@(x)(phideriv(i,x,X,Y)*phideriv(i+1,x,X,Y)),a,b,tol);
 end
 
 A = diag(c1) + diag(c2,1) + diag(c2,-1);
 
 if(N>=4)
     %if you have more then 3 basis functions
     %then you must integrate phi(2i)*phi(2i+2)
     %(for the even basis functions only)
     temp = N-4;
     for(i=1:floor(temp/2)+1)
         a = X(1+i);
         b = X(2+i);
         c3 = quadv(@(x)(phideriv(2*i,x,X,Y)*phideriv(2*i+2,x,X,Y)),a,b,tol);
         A(2*i,2*i+2) = c3;
         A(2*i+2,2*i) = c3;
     end
 end
 A = sparse(A);
 
 %construct F
 F = zeros(N,1);
 
 for(i = 1:N) %n-1
%    high order quadratures to approx (f,phi(i))
     odd = (i+1)/2-1; 
     even = i/2-1;
     if(mod(i,2)~=0) %if odd phi
         a = X(1+odd);
         b = X(2+odd);
     else
         a = X(1+even);
         b = X(3+even);
     end
     F(i)=quadv(@(x)(f(x)*phi(i,x,X,Y)),a,b,tol);
 end
 
 %finally, find U
 U = (A\F);
 U = [0,U',0];
 
 %construct the x-axis, made from weaving together
 %the X and Y partitions:

 bol = 1;
 count = 1;
 for(i=1:N)
     if(bol == 1)
         xAxis(i)=X(count);
         bol = 2;
     else
         xAxis(i)=Y(count);
         bol = 1;
         count = count+1;
     end
 end
 xAxis = [xAxis,Y(length(Y)),1];
 
 hold on;
 plot([xAxis],U,'k');
 plot([0:1/10000:1],sin(pi*[0:1/10000:1]))
 hold off;
 
 err = 0;
 for(i = 1:length(xAxis))
     err = err + (sin(pi*xAxis(i))-U(i))^2;
 end
 err = sqrt(err);
 
 
 
 
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
 %z = sin(pi*x)+exp(2*x);


function z = phi(k,x,X,Y)
%represents the basis function phi
%k is the kth basis function
%x is the location
%X is the X partition (x0,x1,...)
%Y is the Y partition (y1,y2,...)
    
d = (k+1)/2-1;
b = k/2-1;

%only if k is odd!!!
if(mod(k,2)~=0 && x>X(1+d) && x<X(2+d))
    z = -(x-X(1+d))/(Y(1+d)-X(2+d))*((x-X(2+d))/(Y(1+d)-X(2+d)));
elseif(mod(k,2)==0 && x>X(1+b) && x<X(2+b))
    z = (x-X(1+b))/(X(2+b)-X(1+b))*((x-Y(1+b))/(X(2+b)-Y(1+b)));
elseif(mod(k,2)==0 && x>X(2+b) && x<X(3+b))
    z = (x-Y(2+b))/(X(2+b)-Y(2+b))*((x-X(3+b))/(X(2+b)-X(3+b)));
else
    z = 0;
end

function z = phideriv(k,x,X,Y)
%represents the derivative of phi
%k is the kth basis function
%x is the location
%X is the X partition (x0,x1,...)
%Y is the Y partition (y1,y2,...)
    
d = (k+1)/2-1;
b = k/2-1;

%only if k is odd!!!
if(mod(k,2)~=0 && x>X(1+d) && x<X(2+d))
    z = -(2*x-X(2+d)-X(1+d))/((-Y(1+d)+X(1+d))*(Y(1+d)-X(2+d)));
elseif(mod(k,2)==0 && x>X(1+b) && x<X(2+b))
    z = (2*x-Y(1+b)-X(1+b))/((X(1+b)-X(2+b))*(Y(1+b)-X(2+b)));
elseif(mod(k,2)==0 && x>X(2+b) && x<X(3+b))
    z = (2*x-X(3+b)-Y(2+b))/((X(2+b)-Y(2+b))*(X(2+b)-X(3+b)));
else
    z = 0;
end
        


function z = integrand(j,x,part)
z=(f(x).*phi(j,x,part));