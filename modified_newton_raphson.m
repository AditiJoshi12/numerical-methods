%% Aditi Joshi
% 200121022 
% Modified Newton Raphson Method to accomodate multiplicity
% Solving system of non-linear equations using Newton Raphson Method
%% Problem 1
%  The equation $(x − 1.1)^2 *(x + 1) = 0$ has a double root at ξ = 1.1. 
% Write a program to compute an approximate root to ξ by using the standard Newton’s formula as well as
% the modified Newton’s formula. Determine the order of convergence numerically for both the cases. 
% Compute the order of convergence with the formula 
% $p = log10 ( |en+2| / |en+1| )/ log10 ( |en+1| / |en| )$ , n = 0, 1, 2, .
% where en = ξ − xn. Take TOL = 10−3 and x0 = 1.
%% Plotting the function

% x = -2:0.01:3;
% y = (x+1).*(x-1.1).^2;
% 
% plot(x, y)
% title('Plotting (x+1)*(x-1.1)^2')
% xlabel('X')
% ylabel('f(x)')
%% Newton Raphson Method

% defining the function and derivative of function

f = @(x) (x+1)*(x-1.1)^2;
f_prime = @(x) 3*(x-1.1)*(x+0.3);

iter = 0; % initializing iterations
max_iter = 100; % setting limit on iterations
tol = 0.001; % declaring tolerance
x0 = 1; % initial value for root
x1 = x0 - (f(x0)/f_prime(x0)); % better approximated root
x2 = x1 - (f(x1)/f_prime(x1)); % next approximation for calculating order of convergence

order = log10(abs((1.1 - x2)/(1.1 - x1)))/log10(abs((1.1 - x1)/(1.1 - x0))); % order of convergence

fprintf("Newton Raphson Method\n\nIteration\tNewton-Raphson Method\tOrder of Convergence\n")
fprintf('\n%d\t\t\t%5d\t\t\t%5d\n', iter+1, x1, order)

while abs(x1 - x0) > tol % while error is more than tolerance
    if iter < max_iter % if number of iterations less than max iterations
        x0 = x1; % x0 -> x1
        x1 = x0 - (f(x0)/f_prime(x0)); % better approximation
        x2 = x1 - (f(x1)/f_prime(x1)); % next approximation for calculating order of convergence
        order = log10(abs((1.1 - x2)/(1.1 - x1)))/log10(abs((1.1 - x1)/(1.1 - x0))); % order of convergence
        iter = iter + 1; % increase number of iterations
    end
    fprintf('\n%d\t\t\t%5d\t\t\t%5d\n', iter+1, x1, order) % printing results
end

fprintf("\nNumber of Iterations for Newton Raphson Method are %d.\n\n", iter+1)
%% Modified Newton Raphson Method

% The multiplicity of root(m) is 2

f = @(x) (x+1)*(x-1.1)^2;
f_prime = @(x) 3*(x-1.1)*(x+0.3);

iter = 0;
max_iter = 100;
tol = 0.001;
m = 2; % multiplicity of the root
x0 = 1; % initial aproximation of the double root
x1 = x0 - m*(f(x0)/f_prime(x0)); % modified Newton Raphson Method
x2 = x1 - m*(f(x1)/f_prime(x1)); % modified Newton Raphson Method

order = log10(abs((1.1 - x2)/(1.1 - x1)))/log10(abs((1.1 - x1)/(1.1 - x0))); % order of convergence

fprintf("Modified Newton Raphson Method\n\nIteration\tNewton-Raphson Method\tOrder of Convergence\n")
fprintf('\n%d\t\t\t%5d\t\t\t%5d\n', iter+1, x1, order)

while abs(x1 - x0) > tol
    if iter < max_iter
        x0 = x1;
        x1 = x0 - m*(f(x0)/f_prime(x0));
        x2 = x1 - m*(f(x1)/f_prime(x1));
        if x2 ~= 1.1
            order = log10(abs((1.1 - x2)/(1.1 - x1)))/log10(abs((1.1 - x1)/(1.1 - x0)));
        end        
        iter = iter + 1;
    end
    fprintf('\n%d\t\t\t%5d\t\t\t%d\n', iter+1, x1, order)
end

fprintf("\nNumber of Iterations for Modified Newton Raphson Method are %d.\n\n", iter+1)
%% Systems of Non Linear Equations

% Write a MATLAB program to compute an approximate solution of the following nonlinear
% system
% f1(x1, x2) := sin(x1x2) + x1 − x2 = 0,
% f2(x1, x2) := x2 cos(x1x2) + 1 = 0,
% using Newton’s method. Take the starting value [x01, x02] = [1, 2] and use stopping criteria for
% accepting the solution is TOL = 10−3.

%% Plotting the function 

[x, y] = meshgrid(-2:0.1:2);

z1 = sin(x.*y) + x - y;
z2 = y*cos(x.*y) + 1;
surfl(z1)
shading interp

mesh(z2)

%% Solving System of Non Linear Equations

x0 = [1;2]; % initial approximation for root
X1 = [0;0]; % declaring next approximation
iter = 1;
tol = 0.001; % tolerance
h = ones(2, 1); % correction vector

fprintf('Solving Non Linear Systems of Equations\n');
fprintf('\nIteration\tx1\t\t\t\t\tx2\t\t\t\t\tf1(x1,x2)\t\t\tf2(x1,x2)\n');
fprintf('\n%d\t\t\t%d\t\t\t\t\t%d\t\t\t\t\t%d\t\t\t\t\t%d\n\n',iter, x0(1), x0(2), f1(0, 0), f2(0, 0));

while abs(h(1)) + abs(h(2)) > tol
    % defining the Jacobian
    J = [df1_by_dx1(x0(1), x0(2)) df1_by_dx2(x0(1), x0(2)); df2_by_dx1(x0(1), x0(2)) df2_by_dx2(x0(1), x0(2))]; 
    F = [f1(x0(1), x0(2)); f2(x0(1), x0(2))]; % function vector
    h = - J \ F; % correction vector
    X1 = x0 + h; % better approximation
    x0 = X1; % x0 -> x1
    iter = iter+1; % step counter
    
    fprintf('%d\t\t\t%d\t\t%d\t\t%d\t\t%d\n\n',iter, X1(1), X1(2), f1(X1(1), X1(2)), f1(X1(1), X1(2)));
end

% defining functions

function f_1 = f1(x1, x2)
    f_1 = sin(x1*x2) + x1 - x2;
end

function f = f2(x1, x2)
    f = x2 * cos(x1 * x2) + 1;
end

function df1x1 = df1_by_dx1(x1, x2)
    df1x1 = sin(x1*x2) + x1 - x2;
end

function df1x2 = df1_by_dx2(x1, x2)
    df1x2 = x1*cos(x1*x2) - 1;
end

function df2x1 = df2_by_dx1(x1, x2)
    df2x1 = - sin(x1 * x2) * x2^2;
end

function df2x2 = df2_by_dx2(x1, x2)
    df2x2 = cos(x1*x2) - sin(x1 * x2) * x1 * x2;
end