clear all
close all

%%Problem 1
H = [0.0001, 0.001, 0.01, 0.1, 1, 2, 3];
hl = length(H);
for i = 1 : hl
    [y, yy, T] = Lobatto(H(i));
    yE = Euler(H(i));
    figure(i)
    plot (T, y, T, yy, T, yE)
    xlabel('t')
    ylabel('y')
    legend('Lobatto IIID', 'ode45', 'Euler')
    title(['Numerical solution obtained for h = ',num2str(H(i)),' by Lobatto method, function ode45 and Euler method']);
end

%Problem 2
H = logspace(-5, 0, 20);
hl = length(H);

accuracy = zeros(size(H));
accuracyInf = zeros(size(H));
accuracyEuler = zeros(size(H));
accuracyInfEuler = zeros(size(H));

for i = 1 : hl
    h = H(i);
    [y, yy, T] = Lobatto(h);
    accuracy(i) = norm(y - yy) / norm(yy);
    accuracyInf(i) = norm(y - yy, Inf) / norm(yy, Inf);
    
    yEuler = Euler(h);
    accuracyEuler(i) = norm(yEuler - yy) / norm(yy);
    accuracyInfEuler(i) = norm(yEuler - yy, Inf) / norm(yy, Inf);
    
end

figure(8)
loglog(H, accuracy, H, accuracyEuler);
xlabel('h');
ylabel('\bf \delta_{2}');
legend('\bf \delta_{2 Lobatto}', '\bf \delta_{2 Euler}');
title(['The dependence of \delta_{2}(h) on h']);

figure(9)
loglog(H, accuracyInf, H, accuracyInfEuler);
xlabel('h');
ylabel('\bf \delta_{\infty}');
legend('\bf \delta_{\infty Lobatto}', '\bf \delta_{\infty Euler}');
title(['The dependence of \delta_{\infty}(h) on h']);



%Solving using the implicit method Lobatto IID of order 2
function [y, yy, T] = Lobatto(h)
    T = linspace (0, 10, 10/h + 1);
    y = zeros(size(T));
   
    a = 1/2*h;
    A = [1 + a, a;
        -a, 1 + a];
    
    B = zeros(2, 1);
    F = zeros(2, 1);
    
    N = length(y);

    for n = 2 : N
        B = [-y(n-1) + 2 * T(n-1) * exp(-T(n-1) + 2);
            -y(n-1) + 2 * (T(n-1) + h) * exp(-T(n-1) - h + 2)];
        F = A\B;
        y(n) = y(n-1) + h * (1/2*F(1) + 1/2*F(2));
    end
    y = y';
    
    opts = odeset('RelTol',2.22045e-14,'AbsTol',1e-20);
    yy0 = 0;
    [~, yy] = ode45(@(a, x) -x + 2*a*exp(-a+2), T, yy0, opts);
end



%Solving using the explicit Euler method
function [y] = Euler(h)
    T = linspace (0, 10, 10/h + 1);
    y = zeros(size(T));
    N = length(y);

    for n = 2 : N
        y(n) = y(n-1) + h * (-y(n-1) + 2 * T(n-1) * exp(-T(n-1) + 2));
    end
    y = y';
end