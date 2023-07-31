% Apply composite trapezoid and corrected trapezoid rules to the approximation of integral
% integral 0 to 1 x^2*e^−2xdx = 0.0808308960 . . . ,
% and compare your results in light of the expected error theory for both methods, and comment on
% what occurs. Take N = 50.

trap = trapezoidal();
corr_trap = corrected_trapezoidal();

et = 0.0808308960 - trap;
ect = 0.0808308960 - corr_trap;

fprintf('\nIntegral using composite trapezoidal rule is %8f and the error is %d\n', trap, et);
fprintf('\nIntegral using composite corrected trapezoidal rule is %8f and the error is %d\n', corr_trap, ect);


% f(x) = @(x) x*x*exp(-2*x);
% f_prime(x) = @(x) 2*x*exp(x) - 2*x*x*exp(-2*x);

function F = f(x)
    F = x*x*exp(-2*x);
end

function F_der = f_prime(x)
    F_der = 2*x*exp(-2*x) - 2*x*x*exp(-2*x);
end

function T = trapezoidal()
    n = 50;
    h = 1/n;
    T = f(0) + f(1);
    for i = linspace(h,1-h, n-1)
        T = T + 2 * f(i);
    end    
    T = h*T/2;
end

function TC = corrected_trapezoidal()
    n = 50;
    h = 1/n;
    TC = f(0) + f(1);
    for i = linspace(h,1-h, n-1)
        TC = TC + 2 * f(i);
    end    
    TC = h*TC/2;
    TC = TC + (h*h*(f_prime(0)-f_prime(1))/12);
end



