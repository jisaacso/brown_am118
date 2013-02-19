function z = phi(j,x,part)

%part = partition of space (= x in FiniteElet.m)
%x = current spatial location
%j = current phi


if(part(j)<=x && part(j+1)>=x)
    z = (x-part(j))./(part(j+1)-part(j));
elseif(x>part(j+1) && x<part(j+2))
    z = (part(j+2)-x)./(part(j+2)-part(j+1));
else
    z = 0;
end