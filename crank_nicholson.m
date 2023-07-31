%% Problem ii)

% A uniform solid rod of one-half a unit of length is thermally insulated along its length and its initial
% temperature at zero time is 0oC. One end is thermally insulated and the other supplied with heat
% at a steady rate. The temperature at points within the rod satisfy the equation
% ∂U/∂t = ∂^2U/∂x^2 
% for all x ∈ (0, 1/2) and t > 0,
% satisfying the initial condition
% U(x, 0) = 0, 0 ≤ x ≤ 1/2,
% and the boundary conditions
% ∂U/∂x = 0 at x = 0, t > 0,
% ∂U/∂x = 1 at x = 1/2, t > 0.
% Solve the problem numerically using the Crank-Nicolson method with h = 0.1 and r = 1. 
% The analytical solution is given by
% U(x, t) = 2t + 1/2((12x^2 − 1)/6 − 2/π^2 sigma(n: 1 -> ∞) (−1)^n/n^2(e^(−4π^2n^2t))cos(2nπx))
% Print the solution at the time t = 0.01, 0.05, 0.5, 1 for all x = 0, 0.1, 0.2, 0.3, 0.4, 0.5. For a comparison
% with the exact solution U(x, t), take n = 1 in (1).

%% Numerical Solution using Crank-Nicholson Method

h = 0.1;
r = 1;
k = r*h^2;

x = 0:h:0.5;
t = 0:k:1;

n = length(x);
m = length(t);

U = zeros(n, m);

for i = 1:1:m
    U(1, i) = 0;    
end

A = (2 + 2*r)*eye(n-2);

for i = 1:1:n-3
    A(i, i+1) = -r;
    A(i+1, i) = -r;
end

B = (2 - 2*r)*eye(n-2);

for i = 1:1:n-3
    B(i, i+1) = r;
    B(i+1, i) = r;
end

for j = 2:1:m
    U(1, j) = 2*r*U(2, j-1) + (1 - 2*r)*U(1, j-1);
    U(2:n-1,j) = A\B*U(2:n-1,j-1);
    U(n, j) = 2*h*r + 2*r*U(n-1, j-1) + (1 - 2*r)*U(n, j-1);
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
% colormap 'summer'
title('Numerical Solution Using Classical Explicit Scheme');
xlabel('X');
ylabel('T');
zlabel('U(x, t)');

%% Printing solutions

fprintf('\nX\t\tT\t\tExact\t\tExplicit\n');

for x = [0, 0.1, 0.2, 0.3, 0.4, 0.5]
    for t = [0.01, 0.05, 0.5, 1]
        i = int16((x/h)) + 1;
        j = int16((t/k)) + 1;
        exact = U_exact(x, t);
        numer = U(i, j);
        fprintf('\n%0.2f\t%0.2f\t%0.6f\t%0.6f\n',x, t, exact, numer);        
    end
end

%% Defining Exact Solution

function u = U_exact(x, t)
    u = 2*t + x.^2 + (-1/12) + (1/pi^2)*exp(-4*pi^2*t).*cos(2*pi*x);
end