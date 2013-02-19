%%%%%%%%%%%%%%%%%
%@Joseph Isaacson%
%%%%%%%%%%%%%%%%%%

function hw4_3(t0,T,u0,k,f)

% hw4_3 uses Forward and Backward Euler method
% To solve a 2x2 system of ODEs
% whose Right Hand Side (RHS) is stored
% within the .m file "hw4_3ODE"
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
U = zeros(2,N);
U(1,1)= u0(1);
U(2,1)= u0(2);
u0 = u0';

for(n=1:N-1)
   
    %forward Euler:
    %U(:,n+1) = U(:,n)+k*f(U(:,n),t(n));
    
    %backward Euler:
    U(:,n+1) = fsolve(@(utemp) utemp-U(:,n)-k*f(utemp,t(n+1)), u0);
    %U(:,n+1)
end
%print error
exact = [exp(-1);exp(-1000)];
abs(norm(U(:,N))-norm(exact))
