function hw2_1(a,b,k,f)

% hw2_1 uses Simpson's method for approximating
% a definite integral whose integrand is stored
% within the .m file "hw2_1integrand"
% INPUTS:
%    a  = lower bound on integral
%    b  = upper bound on integral
%    k  = step size
%    f  = integrand (from command line,
%         MUST be entered as: @hw2_1integrand)
% OUTPUT:
%    U = solution of integral


% Initialization:
N=(b-a)/k;              %N = num of steps
U = 0;                  %approx of integral
t = linspace(a,b,N);    %pre-allociate space for t
                        %to make Simpson's formula
                        %more elegant looking (note:
                        %if memory is an issue, just
                        %update t every time step)

% Main Code:
for(n=1:N-1)
    U = U+(t(n+1)-t(n))/6*(f(t(n))+4*f((t(n)+t(n+1))/2)+f(t(n+1)));
end
% Print final integral value:
U

