%% Problem 1
% Consider the one-dimensional heat equation $\frac{\partial U}{\partial t} 
% = \frac{\partial^2 U}{\partial x^2}$ for (x, t) ∈ (0, 1) × (0, 0.1]. Use classical 
% explicit scheme to solve the above equation, where the initial and boundary 
% data are taken from the exact solution U(x, t) = exp(−π 2 t) sin πx. Perform 
% the following experiments: (i) Plot the numerical solutions with h = 0.1, k 
% = .005 against the exact solution; (ii) Study the convergence of numerical solutions 
% when k/h^2 > 1/2.

%% Part 1 

x = 0:0.1:1;
t = 0:0.005:0.1;

h = 0.1;
k = 0.005;

n = length(x);
m = length(t);

r = k/h^2;

U = zeros(n, m);

for j = 1:1:m
    U(1, j) = U_exact(0, k*j - k);
    U(n, j) = U_exact(1, k*j - k);
end

for i = 1:1:n
    U(i, 1) = U_exact(i*h - h, 0);    
end

for j = 1:1:m-1   
    for i = 2:1:n-1 
        U(i, j+1) = r*U(i+1, j) + (1-2*r)*U(i, j) + r*U(i-1, j);
    end
end

[X, T] = meshgrid(x, t);

surf(x, t, U_exact(X, T))
shading flat;
title('Exact Solution');
xlabel('X');
ylabel('T');
zlabel('U(x, t)');

figure();

s = mesh(X, T, U');
s.FaceColor = 'flat';
title('Numerical Solution Using Classical Explicit Scheme');
xlabel('X');
ylabel('T');
zlabel('U(x, t)');

%% Part 2;

x = 0:0.01:1;
t = 0:0.005:0.1;

h = 0.01;
k = 0.005;

n = length(x);
m = length(t);

r = k/h^2;

U = zeros(n, m);

for j = 1:1:m
    U(1, j) = U_exact(0, k*j - k);
    U(n, j) = U_exact(1, k*j - k);
end

for i = 1:1:n
    U(i, 1) = U_exact(i*h - h, 0);    
end

for j = 1:1:m-1   
    for i = 2:1:n-1 
        U(i, j+1) = r*U(i+1, j) + (1-2*r)*U(i, j) + r*U(i-1, j);
    end
end

[X, T] = meshgrid(x, t);

figure();

s = mesh(X, T, U');
s.FaceColor = 'flat';
title('Numerical Solution Using Classical Explicit Scheme (but with r > 0.5)');
xlabel('X');
ylabel('T');
zlabel('U(x, t)');

%% Defining Exact Solution

function u = U_exact(x, t)
    u = exp(-pi^2*t).*sin(pi*x);
end