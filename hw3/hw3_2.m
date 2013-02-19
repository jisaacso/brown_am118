%%%%%%%%%%%%%%%%%
%@Joseph Isaacson%
%%%%%%%%%%%%%%%%%%

function hw3_2(t0,T,u0,k,f)

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
%N=(T-t0)/k;
N = 1/k;
t = linspace(t0,T,N);
U = zeros(1,N);
U(1)= u0;
U(2) = U(1) + k*f(U(1),t(1));

for(n=1:N-2)
    
    U(n+2) = 3*U(n+1)-2*U(n)-k*f(U(n),t(n));

end
U(N)-exp(-1)