function z = rho(k,x,X,Y)

d = (k+1)/2-1;


% x0 = 0;
% y1 = .5;
% x1 = 1;

%only if k is odd!!!
if(mod(k,2)~=0 && x>X(1+d) && x<X(2+d))
    %z = (x-x0)/(y1-x0)*((x-x1)/(y1-x1))
    z = -(x-X(1+d))/(Y(2+d)-X(2+d))*((x-X(2+d))/(Y(2+d)-X(2+d)));
else
    z = 0;
end