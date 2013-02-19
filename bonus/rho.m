function z = rho(k,x,X,Y)

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