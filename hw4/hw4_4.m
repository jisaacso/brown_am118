%%%%%%%%%%%%%%%%%
%@Joseph Isaacson%
%%%%%%%%%%%%%%%%%%

function hw4_4(t0,T,u0,k,f,d,kbar)

% hw4_4 uses euler's method with adaptive
% time step sizes to solve an ODE
% whose Right Hand Side (RHS) is stored
% within the .m file "hw4_4ODE"
% INPUTS:
%    t0   = initial time
%    T    = final time
%    u0   = initial value of ODE
%    k    = step size
%    f    = RHS of the ODE (from command line,
%           MUST be entered as: @hw1_1ODE)
%    d    = Tolerance
%    kbar = Smallest Step Size allowed
% OUTPUT:
%    U = solution of ODE at time T


% initialize U, t, and define N = num time steps
N = 1/k;
U = u0;
t = t0;

while(t < N)
    Wstar = U + k/2*f(U,t);           %2nd order runga-kutta
    W = U + k*f(Wstar,t+k/2);         %method acts as exact ans
    
    V = U+k*f(U,t);                   %euler method to look at
                                      %estimated soln to ode
                                          
    errorNextStep = abs(W-V);         %find error between 2 methods
    
    %if the error > tolerance
    %AND the current step size is not too small
    if(errorNextStep > k*d && (k/2) > kbar) 
        k = k/2;                      %decrease step size
    else
        %else, if error is too small (program too slow)
        if(errorNextStep <= k*d/10)
            k = 2*k;                  %increase step size
        end
        U = W;                        %set current approx to RK2's output
    end
    t = t+k;                          %increment time

end
error = abs(U-exp(-1000))