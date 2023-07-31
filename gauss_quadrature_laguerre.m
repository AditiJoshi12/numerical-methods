% Consider two-point Gauss quadrature rule (n = 2)
% Z ∞,0 e^−xf(x)dx ≈ W1f(x1) + W2f(x2).
% Note that, here the weight function is w(x) = e^−x, so the orthogonal polynomial family that one needs
% to use is the Laguerre family. The Gauss points x1 and x2 are the roots of the second-order Laguerre
% polynomial L2(x) = 1/2(x^2 - 4x + 2). The roots(the Gauss points) are x1 = 2 - sqrt(2) = 0.5857864376, 
% x2 = 2 + sqrt(2) = 3.414213562. The weights are W1 = 0.8535533903, W2 = 0.1464466092 (Can you
% verify these two values?). Compute the approximate value of the integral for f(x) = x^3.

%% 

x1 = 0.5857864376;
x2 = 3.414213562;
W1 = 0.8535533903;
W2 = 0.1464466092;

I = W1*(x1^3) + W2*(x2^3);

fprintf('\nThe approximation of the given integral using Gauss Quadrature Rule is = %0.8f\n', I);
