function hw2_4order(a,b,f)

% hw2_4 finds the order of the 2 and 3 point Gauss quad. methods 
% for approximating a definite integral whose 
% integrand is stored within the .m file "hw2_1integrand"
% INPUTS:
%    a  = lower bound on integral
%    b  = upper bound on integral
%    k  = step size
%    f  = integrand (from command line,
%         MUST be entered as: @hw2_1integrand)
% OUTPUT:
%    U = solution of integral

numExps = 10;
e = zeros(1,numExps);
for(exponent=1:numExps)
    k = 1/2^exponent;


    % Initialization:
    N  =(b-a)/k;                   %N = num of steps
    U  = zeros(1,N);               %approx of integral

    b2 = [1,1];                    %2-point gauss weights
    x2 = [-1/sqrt(3),1/sqrt(3)];   %2-point gauss points
    b3 = [8/9,5/9,5/9];            %3-point gauss weights
    x3 = [0,sqrt(3/5),-sqrt(3/5)]; %3-point gauss points

    t  = linspace(a,b,N);          %pre-allociate space for t
                                   %to make gauss quad formula
                                   %more elegant looking (note:
                                   %if memory is an issue, just
                                   %update t every time step)
              

    % Main Code:
    for(n=1:N-1)
    % 2 point Gauss Quad:
        for(i=1:2)
            U = U+b2(i)*f(k*x2(i)/2+(t(n)+t(n+1))/2);
        end
    
        % 3-point Gauss Quad:
%         for(i=1:3)
%             U(n)= U(n)+b3(i)*f(k*x3(i)/2+(t(n)+t(n+1))/2);
%         end
%         U(n+1)=U(n);
   
    end
    U = U*k/2;
    
    %find error of Gauss for given step size k:
    e(exponent) = abs(U(N)-(1-exp(-1)));
end

% Calculate Order of Method:
for(i=1:numExps-1)
    order = log(e(i)/e(i+1))/log(2)
end
