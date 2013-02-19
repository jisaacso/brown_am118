%%%%%%%%%%%%%%%%%
%@Joseph Isaacson%
%%%%%%%%%%%%%%%%%%

function hw1_1(t0,T,u0,k,f)

% hw1_1 uses varying methods to solve an ODE
% whose Right Hand Side (RHS) is stored
% within the .m file "hw1_1ODE"
% INPUTS:
%    t0 = initial time
%    T  = final time
%    u0 = initial value of ODE
%    k  = step size
%    f  = RHS of the ODE (from command line,
%         MUST be entered as: @hw1_1ODE)
% OUTPUT:
%    U = solution of ODE at time T


%Uncomment below to find order of methods:
for(i=1:6)
k = 1/(2^i);

% initialize U, t, and define N = num time steps
N=(T-t0)/k;
U = zeros(1,N);

U(1) = u0;
t = linspace(t0,T,N);

for(n=1:N-1)
    %Euler's Method:
    %U(n+1)=U(n)+k*f(U(n),t(n));
    
    %Backward Euler's:
    %U(n+1) = fsolve(@(utemp) utemp-U(n)-k*f(utemp,t(n+1)), u0);
    
    %Trapezoidal:
    %U(n+1) = fsolve(@(utemp) utemp-U(n)-k/2*(f(U(n),t(n))+f(utemp,t(n+1))), u0);
    
    %Midpoint:
    U(2) = U(1)+k*f(U(1),t(1));
    U(n+2) = U(n)+2*k*f(U(n+1),t(n+1));
     
    %2nd Order Runge-Kutta
%     Ustar  = U(n) + k/2*f(U(n),t(n));
%     U(n+1) = U(n) + k*f(Ustar,t(n)+k/2);
    
    %4th Order Runge-Kutta
%     Y1 = U(n);
%     Y2 = U(n) + k/2*f(Y1,t(n));
%     Y3 = U(n) + k/2*f(Y2,t(n)+k/2);
%     Y4 = U(n) + k/2*f(Y3,t(n)+k/2);
%     U(n+1) = U(n) + k/6*(f(Y1,t(n))+2*f(Y2,t(n)+k/2)+2*f(Y3,t(n)+k/2)+f(Y4,t(n)+k));

end

%Calculate and Plot error
%(comment all code below when
%finding the order of a method)
% e(1)=0;
% for(i=2:N)
%     e(i) = abs(U(i)-exp((i)*k)+.5*exp(-(i)*k));
% end
% plot(t,e,'k',t,e1,'--k');
% legend('4th order error', '2nd order error');
% title('R-K 4th order error plot');

%calculate order of method:
%(comment when plotting error above)
e(i) = abs(U(i)-exp((i)*k)+.5*exp(-(i)*k));
end
p = log(e(5)/e(6))/log(2)







