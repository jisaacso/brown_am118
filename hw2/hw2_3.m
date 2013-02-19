%%%%%%%%%%%%%%%%%
%@Joseph Isaacson%
%%%%%%%%%%%%%%%%%%

function hw2_3(t0,T,u0,k,f)

% hw2_3 uses varying methods to solve an ODE
% whose Right Hand Side (RHS) is stored
% within the .m file "hw2_3ODE"
% INPUTS:
%    t0 = initial time
%    T  = final time
%    u0 = initial value of ODE
%    k  = step size
%    f  = RHS of the ODE (from command line,
%         MUST be entered as: @hw1_1ODE)
% OUTPUT:
%    U = solution of ODE at time T


% initialize U, t, and define N = num time steps
N=(T-t0)/k;
U = u0;
t = linspace(t0,T,N);

for(n=1:N-1)
    
    %Classical 3rd Order Runge-Kutta
    Y1 = U;
    Y2 = U + k/2*f(Y1,t(n));
    Y3 = U - k*f(Y1,t(n))+2*k*f(Y2,t(n)+k/2);
    U = U + k/6*(f(Y1,t(n))+4*f(Y2,t(n)+k/2)+f(Y3,t(n)+k/2));

end
U