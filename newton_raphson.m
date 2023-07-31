%% |*Aditi Joshi*|
% |*200121022*|
% 
% |*MA311M Scientific Computing*|
%% |_*Fixed Point Method, Bisection Method, Newton Raphson Method*_|
%% Problem 1
% Write a MATLAB program implementing the Bisection method and Newton-Raphson 
% method (NRM) to compute the smallest positive root of the equation $f(x) = 0$, 
% where $f(x) = e ^ {-x} - sin x.$
%% 
% * Use the Bisection method to compute the approximate root of $f(x) = 0$ correct 
% to five decimal places. Note the number of iterations is needed to achieve this 
% accuracy.
% * Print the approximate solutions. 
%% Plotting the function

x = 0:0.001:1 ;
y = exp(-x) - sin(x) ;

plot(x, y)
title('Plotting e^{-x} - sin(x)')
xlabel('X')
ylabel('f(x)')
%% Bisection Method

f = @(x) 1/(3*x - 1); %defining the function

num_iter = 0;
a = 0;
b = 1;
tol = 0.001;
c = (a + b)/2; % initializing midpoint

fprintf("\nIteration\tBisection Method\t\n");
fprintf('\n%d\t\t\t\t%5d\n', num_iter, c);

while abs(b-a)>tol
    if f(a)*f(c) < 0
        b = c;
    else
        a = c;
    end
    num_iter = num_iter + 1;
    c = (a + b)/2;
    fprintf('\n%d\t\t\t\t%5d\n', num_iter, c);
end

fprintf('Number of Iterations for Bisection Method are %d.\n', num_iter);
%% Newton Raphson Method

f_prime = @(x) -exp(-x) - cos(x);

a = 0;
b = 1;
iter = 0;
max_iter = 1000000;
tol = 0.00001;
x0 = 0.6;
x1 = x0 - (f(x0)/f_prime(x0));

fprintf("Newton Raphson Method\nIteration\tNewton-Raphson Method\t\n")

while abs(x1 - x0) > tol
    if iter < max_iter
        x0 = x1;
        x1 = x0 - (f(x0)/f_prime(x0));
        iter = iter + 1;
    end
    fprintf('\n%d\t\t\t%5d\n', iter, x1)
end

fprintf("Number of Iterations for Newton Raphson Method are %d.\n", iter)
%% Problem 2
% Write a program which uses *fixed-point iteration* to find the zero of the 
% function $f(x) = x^2 - x - 2$ in the interval $[0, 7]$. 
% 
% Take the iteration function $g(x) = \surd  (x + 2)$ and the starting point 
% $x_0 = 0$. 
% 
% Terminate the programme when the following error tests are satisfied $\frac{|x_n 
% - x_{n-1}|} {|x_n|} < tol$ or $|f(x)| < tol$,
%% Plotting the Function
% 
% X = 0:0.01:7;
% Y = X.^2 - X - 2;
% 
% plot(X, Y)
% title('Plotting x^{2} - x - 2')
% xlabel('X')
% ylabel('f(x)')
%% Fixed Point Method
% $x_{n+1} = g(x_{n})$

f2 = @(x) x*x - x - 2;
g = @(x) sqrt(x+2); %defining the function

x0 = 0;
tol = 0.001;
x1 = g(x0);
iter = 0;
err = abs(2-x1);
fprintf("Fixed Point Iteration Method\nIteration\t\tx_n\t\tf(x_n)\t\t\tError\n")
fprintf("x%d\t\t%5d\t%5d\t%5d\n", iter+1, x1, f2(x1), err);

while abs(x1-x0)/x1 > tol
    x0 = x1;
    x1 = g(x0);
    err = abs(2 - x1);
    iter = iter + 1;    
    fprintf("x%d\t\t%5d\t%5d\t%5d\n", iter+1, x1, f2(x1), err);
end

fprintf("\n")