function AB(t0,T,u0,k,f)

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AB uses the 2-step Adam-Bashforth method     %
% to solve any ODE of the form u(t)'=f(u(t),t) %
% whose Right Hand Side, f(u(t),t) is stored   %
% within the .m file "AB_ODE"                  %
% INPUTS:                                      %
%    t0 = initial time                         %
%    T  = final time                           %
%    u0 = initial value of ODE                 %
%    k  = step size                            %
%    f  = RHS of the ODE (from command line,   %
%         MUST be entered as: @MIDTERM_ODE1)   %
%                                              %
%    sample command line call:                 %
%    AB(0,1,1,1/100,@MIDTERM_ODE1);            %
%                                              %
% OUTPUT:                                      %
%    U = solution of ODE at time T             %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initialization::
N = (T-t0)/k;         %N = number of time steps
U = zeros(1,N);       %U is approx ODE soln. at every time step
U(1) = u0;            %U0 = u0; (initially given)
t = linspace(t0,T,N); %t is a vector of the times we loop over


%We now use Euler Approximation to find u(t0+k):
U(2) = U(1)+k*f(U(1),t(1));
for(n=1:N-2)
    %2 Step Adam Bashforth Method:
    U(n+2)=U(n+1)+k/2*(-f(U(n),t(n))+3*f(U(n+1),t(n+1)));
end
U