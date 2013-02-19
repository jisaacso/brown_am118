function [error, p]=RK_4_Error(M)

% error will be a vector containing the errors
% p will be a vector containing the approximate order
% M is how many k's we will take

for i=1:M
    k=1/2^i;  %the step size
    appsol=RK_4(0,1,1,k,@MIDTERM_ODE1); % this is the euler solution for k
    error(i)=abs(appsol-exp(-1)); % computes the error
    
    % in this problem the exact solution is exp(1); however, the exact
    % solution for your homework is not exactly this
    
end

for i=1:M-1
   p(i)=log(error(i)/error(i+1))/log(2); % approximate order
end