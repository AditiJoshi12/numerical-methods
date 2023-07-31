% Use the Gaussian quadrature rule
% Integral -1 to 1 f(x)dx ≈ 5/9f(−sqrt(3/5)) + 8/9f(0) + 5/9f(sqrt(3/5)
% to evaluate the integral Integral 0 to 4 (sin t)/t dt.
%% 

% 0 -> 4 --> -1 -> 1
% x -> (t-2)/2


I = (5*f(-sqrt(3/5)) + 8*f(0) + 5*f(sqrt(3/5)))/9;

fprintf('\nThe approximation of the given integral using Gauss Quadrature Rule is = %0.8f\n', I);

function fun = f(t)
    fun = sin(2*t + 2)/(t + 1);
end

