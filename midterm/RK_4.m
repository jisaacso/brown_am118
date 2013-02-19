function Output = RK_4(t0,T,u0,k,f)

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AB uses the 4-stage Runga-Kutta method       %
% to solve any ODE of the form u(t)'=f(u(t),t) %
% whose Right Hand Side, f(u(t),t) is stored   %
% within the .m file "RK4_ODE"                 %
% INPUTS:                                      %
%    t0 = initial time                         %
%    T  = final time                           %
%    u0 = initial value of ODE                 %
%    k  = step size                            %
%    f  = RHS of the ODE (from command line,   %
%         MUST be entered as: @MIDTERM_ODE1)   %
%                                              %
%    sample command line call:                 %
%    AB(0,1,1,1/100,@MIDTERM1_ODE);            %
%                                              %
% OUTPUT:                                      %
%    U = solution of ODE at time T             %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initialization::
N = (T-t0)/k;         %N = number of time steps
U = zeros(1,N);       %U is approx ODE soln. at every time step
U(1) = u0;            %U0 = u0; (initially given)
t = linspace(t0,T,N); %t is a vector of the times we loop over


for(n=1:N-1)
    %4-stage Runga-Kutta Method:
    Y1 = U(n);
    Y2 = U(n) + k/2*f(Y1,t(n));
    Y3 = U(n) + k/2*f(Y2,t(n)+k/2);
    Y4 = U(n) + k/2*f(Y3,t(n)+k/2);
    U(n+1) = U(n)+k/6*(f(Y1,t(n))+2*f(Y2,t(n)+k/2)+...
                2*f(Y3,t(n)+k/2)+f(Y4,t(n)+k));
end
Output = U(length(U));