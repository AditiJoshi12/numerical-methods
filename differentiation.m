% Approximate the derivative of f(x) = x log(1 +x)^2 at x = 1 numerically by using forward, backward,
% and central difference formulae i.e.,
% f'(x) ≈ (f(x + h) − f(x))/h,
% f'(x) ≈ (f(x) − f(x − h))/h,
% f'(x) ≈ (f(x + h) − f(x − h))/2h
% with the step size h = 0.1, 0.01 and 0.001. Compute the error between the numerical approximations
% and exact derivatives. Plot the computational error versus the step size.
%% 

fprintf('\nStep Size\t\tApproximate\t\tExact\t\t\tError\n\n')

H = [0.1, 0.01, 0.001];

f = @(x) x*log(1+x^2);
f_prime = @(x) log(1 + x^2) + (2*x^2/(1+x^2));

errfds = zeros(3);
x = 1:3;

fprintf('Forward Difference Scheme\n\n')

for i = 1:3
    h = H(i);
    fder_approx = (f(1+h) - f(1))/h;
    fder_exact = f_prime(1);
    error = fder_exact - fder_approx;
    errfds(i) = abs(error);
    fprintf('%d\t%d\t%d\t%d\n', h, fder_approx, fder_exact, error);
end

semilogy(x, (errfds), '+:');

hold on;

fprintf('\nBackward Difference Scheme\n\n')

errbds = zeros(3);

for i = 1:3
    h = H(i);
    fder_approx = (f(1) - f(1-h))/h;
    fder_exact = f_prime(1);
    error = fder_exact - fder_approx;
    errbds(i) = abs(error);
    fprintf('%d\t%d\t%d\t%d\n', h, fder_approx, fder_exact, error);
end

semilogy(x, (errbds), 'o-.');

fprintf('\nCentral Difference Scheme\n\n')

errcds = zeros(3);

for i = 1:3
    h = H(i);
    fder_approx = (f(1+h) - f(1-h))/(2*h);
    fder_exact = f_prime(1);
    error = fder_exact - fder_approx;
    errcds(i) = abs(error);
    fprintf('%d\t%d\t%d\t%d\n', h, fder_approx, fder_exact, error);
end

semilogy(x, (errcds), '*--');

hold off;

legend('FDS', 'BDS', 'CDS', 'Location', 'northeast')