x0 = 0;
y0 = [1; 1; 1];

N = [10, 20, 40, 80];

lambda = [-1, 0, 1];
fprintf('Eulers Method and Runge Kutta Method for lambda : %d  %d  %d\n\n', lambda(1), lambda(2), lambda(3));
%fprintf('N\t\tEM\t\tRK4\t\tErr(EM)\t\tErr(RK4)\n\n');

for n = N 
    h = 1/n;
    fprintf('N: %d\n', n);
    exx = exact(1, lambda);
    fprintf('EXACT: %d  %d  %d\n', exx(1), exx(2), exx(3));
    eu = eulers(n, lambda);
    fprintf('EU:    %d  %d  %d\n', eu(1), eu(2), eu(3));
    rk4 = RK4(n, lambda);
    fprintf('RK4:   %d  %d  %d\n', rk4(1), rk4(2), rk4(3));
    err_euf = abs(exx - eu);
    err_eu = abs(err_euf(1) + err_euf(2) + err_euf(3));
    err_rk4f = abs(exx - rk4);
    err_rk4 = abs(err_rk4f(1) + err_rk4f(2) + err_rk4f(3));
    %fprintf('\nN: %d EU: %d\nRK4 : %d\nERROR(EU) : %d\nERROR(RK4) : %d\n', n, eu, rk4, err_eu, err_rk4);    
end

function F = f(x, y, l)
    A = 0.5*[[l(2) + l(3), l(3) - l(1), l(2) - l(1)]; [l(3)-l(2), l(1) + l(3), l(1) - l(2)]; 
        [l(2) - l(3), l(1) - l(3), l(2) + l(1)]];
    F = A*y;
end

function eu = eulers(n, l)
    y1 = zeros(n);
    y2 = zeros(n);
    y3 = zeros(n);
    y1(1) = 1;
    y2(1) = 1;
    y3(1) = 1;
    h = 1/n;
    for i = 1:1:n-1
        xi = i*h;
        y = [y1(i);y2(i);y3(i)];
        Ay = f(xi, y, l);
        y1(i+1) = y1(i) + h*Ay(1);
        y2(i+1) = y2(i) + h*Ay(2);
        y3(i+1) = y3(i) + h*Ay(3);
    end
    eu = [y1(n);y2(n);y3(n)];
end

function Rk4 = RK4(n, l)
    y1 = zeros(n);
    y2 = zeros(n);
    y3 = zeros(n);
    y1(1) = 1;
    y2(1) = 1;
    y3(1) = 1;
    h = 1/n;
    for i = 1:1:n-1
        xi = i*h;
        y = [y1(i);y2(i);y3(i)];
        k1 = f(xi, y, l);
        k2 = f(xi + 0.5*h, y + k1*0.5*h, l);
        k3 = f(xi + 0.5*h, y + k2*0.5*h, l);
        k4 = f(xi + h, y + k3*h, l);
        Y = y + (h*(k1+2*k2+2*k3+k4)/6);        
    end
    Rk4 = Y;
end

function ex = exact(x, l)
    ex(1) = -exp(l(1)*x) + exp(l(2)*x) + exp(l(3)*x);
    ex(2) = +exp(l(1)*x) - exp(l(2)*x) + exp(l(3)*x);
    ex(3) = +exp(l(1)*x) + exp(l(2)*x) - exp(l(3)*x);
end