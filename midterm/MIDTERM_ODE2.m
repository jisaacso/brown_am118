function [z] = MIDTERM_ODE2(u,t)
if(0<=t<=pi)
    z=cos(t);
elseif(pi<t<=5)
    z=sin(t);
else
    output = 'ERROR, TIME MUST RANGE BETWEEN: 0<=t<=5'
end
