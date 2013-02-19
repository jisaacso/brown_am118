function [error, p]=approxorder(M)

% error will be a vector containing the errors
% p will be a vector containing the approximate order
% M is how many k's we will take

for i=1:M
    k=1/2^i;  %the step size
    error(i)=FiniteElet(0,1,k); % this is the FE error for k
end

for i=1:M-1
   p(i)=log(error(i)/error(i+1))/log(2); % approximate order
end