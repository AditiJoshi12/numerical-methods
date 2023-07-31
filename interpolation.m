%% Aditi Joshi - 200121022
% Lab 03
% Interpolation

%% Quadratic Interpolation

%  The error function is an important special function in applied mathematics, with applications to
% probability theory and the solution of heat conduction problems. The formal definition of the error
% function is as follows.
% erf(x) = 2/√π Int(0 to x)(e^-t^2)dt
% Construct the quadratic interpolating polynomial to the error function using the above data at the
% nodes x0 = 0, x1 = 0.5, and x2 = 1.0. Plot the polynomial and the data in the table and comment
% on the observed accuracy.

%% Plotting the function

x = 0:0.01:1;
y = erf(x);

figure()
plot(x, y);
title('erf(x)')
xlabel('X')
ylabel('erf(x)')

%% Forming the quadratic interpolating polynomial

% f(x) = @(x) - 4*0.52049987781305*x*(x-1) + 2*0.84270079294971*x*(x-0.5);

x = 0:0.01:1;
y1 = f(x);
figure()
plot(x, y1);
hold on
plot(x, y, 'Color', 'r');
plot(0, erf(0), '-o', 'Color', '#EDB120')
plot(0.5, erf(0.5), '-o', 'Color', '#EDB120');
plot(1, erf(1), '-o', 'Color', '#EDB120');
hold off

legend('Quadratic Interpolating Polynomial', 'erf(x)', 'Nodes', 'Location', 'southeast');
title('Quadratic Interpolating Polynomial')

% We observe that even with polynomial of degree at most two, accuracy is very high.

%% Newton's Interpolation

% Let f(x) = e^x . Define pn(x) to be the Newton interpolating polynomial for f(x), using n + 1 equally
% spaced nodes on the interval [−1, 1]. Thus, we are taking higher and higher degree polynomial approximations 
% to the exponential function. Write a program that computes pn(x) for n = 2, 4, 8, 16, 32, and
% which samples the error f(x) − pn(x) at 501 equally spaced points on [−1, 1]. Record the maximum
% error as found by the sampling, as a function of n, i.e., define En as
% En = max   |f(tk) − pn(tk)|
%    0≤k≤500
% where tk = −1 + 2k/500 . 
% Print E2, E4, E8, E16, E32 and plot En versus n.

%% Newton's Interpolating Polynomial

X_Newton = -1:0.01:1;

fprintf('\nNewtons Interpolating Polynomial\n\nN\t\tE_N\n');

figure()
plot(X_Newton, exp(X_Newton));
hold on
for n = [2, 4, 8, 16, 32]
    plot(X_Newton, newton(X_Newton, n))
    fprintf('\n%d\t\t%d\n', n, max_error(n))
end   
hold off

legend('Actual', '2', '4', '8', '16', '32', 'Location', 'northwest');
title('Newtons Interpolating Polynomials')

figure()
n = [2,4,8,16,32];
error_f = [max_error(2), max_error(4), max_error(8), max_error(16), max_error(32)];
plot(n, error_f, '-o')
title('Error in Newtons Interpolating Polynomial')
xlabel('n')
ylabel('Error')

%% Lagrange's Interpolating Polynomial

% Write a program to approximate the function
% f(x) = 1/1 + x^2, −5 ≤ x ≤ 5, (Runge’s example),
% using the points xi = −5 + 10*i/8 , i = 0, 1, 2, . . . , 8 by Lagrange’s interpolating polynomial. Plot the
% polynomial against the exact function

X_Lagrangian = -5:0.01:5;

y_Lagrangian = lagrangian(X_Lagrangian);

figure()
plot(X_Lagrangian, y_Lagrangian, 'Color', 'r')
hold on
y_q3 = Q3(X_Lagrangian);
plot(X_Lagrangian, y_q3, 'Color', 'b')
hold off

legend('Lagranges Interpolating Polynomial', 'exp(x)', 'Location', 'southeast');
title('Lagranges Interpolating Polynomial')

%% Helper Functions

function y = f(x)
    y = - 4*0.52049987781305*x.*(x-1) + 2*0.84270079294971*x.*(x-0.5);
end

function f = newton(x, n)
    h = 2/n;
    x0 = -1:h:1;
    
    divided_difference = div_diff(x0, exp(x0));
    
    f = 0;    
    for i = 0:n
        s = divided_difference(1, i+1);
        
        for j = 0:i
            xj = -1 + h*j;
            s = s.*(x - xj);
        end
        f = f + s;
    end    
        
end

function errmax = max_error(n)
    errmax = 0;
    for k = 0:500
        tk = -1 + (k/250);
        err = abs(exp(tk) - newton(tk, n));
        if err > errmax
            errmax = err;
        end
    end
end

function mat = div_diff(x, y)
    m = length(x);
    
    mat = zeros(m, m);
    
    mat(:,1) = y(:); 
    for j = 2:m
        for i = 1:m-j+1
            mat(i,j) = (mat(i+1,j-1)-mat(i,j-1))/(x(i+j-1)-x(i));
        end
    end
    
end

function f_l = lagrangian(x)
    f_l = 0;
    for i = 0:1:8
        xi = -5 + (10*i/8);
        fs_l = Q3(xi);
        for j = 0:1:8
            xj = -5 + (10*j/8);
            if i ~= j
                fs_l = fs_l.*(x - xj)./(xi - xj);
            end
        end
        f_l = f_l + fs_l;
    end
end

function q3 = Q3(x)
    q3 = 1./(1+x.^2);
end
