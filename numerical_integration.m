% The function g(x) is defined by
% g(x) = integral 0 to x e^−x^2dx.
% Write a program for composite rectangle rule (RC), trapezoidal rule (TC) and Simpson’s rule (SC)
% to evaluate g(1) with N = 50, 100, 200 subdivisions. Compare the results with the correct value
% g(1) = 0.74682413 and print the approximate values for RC, TC, SC and the corresponding errors
% ER, ET , ES as per the format shown below.

N = [50, 100, 200];

fprintf('\nN\tRC\t\t\tTC\t\t\tSC\t\t\tER\t\t\tET\t\t\tES\n\n');


for n = N
    rc = rectangle(n);
    er = abs(0.74682413 - rc);
    tc = trapezoidal(n);
    et = abs(0.74682413 - tc);
    sc = simpsons(n);
    es = abs(0.74682413 - sc);
    fprintf('\n%d\t%8f\t%8f\t%8f\t%8f\t%8f\t%8f\n\n', n, rc, tc, sc, er, et, es);
end


function RC = rectangle(n)
    h = 1/n;
    RC = 0;
    for i = linspace(0,1, n+1)
        RC = RC + h * exp(-i*i);
    end    
end

function TC = trapezoidal(n)
    h = 1/n;
    TC = h + h * exp(-1);
    for i = linspace(h,1-h, n-1)
        TC = TC + 2 * h * exp(-i*i);
    end    
    TC = TC/2;
end

function SC = simpsons(n)
    h = 1/n;
    SC = 0;
    
    for i = 1:1:n/2
        a = h*(2*i - 2);
        mp = h*(2*i-1);
        b = h*(2*i);
        SC = SC + exp(-a*a) + 4*exp(-mp*mp) + exp(-b*b);
    end
    
    SC = SC*h/3;
end

















