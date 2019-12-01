clear all
close all
clc

f = @(x) -sin(pi*x).*exp(-x);
Ns = [5, 10, 15];
Ks = [2, 5, 8];
x1 = linspace(-1, 1, 100);

%%Problem 1
for i = 1 : 3
    N = Ns(i);
    y = f(x1);
    
    figure (i)
    plot(x1, y);
    hold on
    grid on
    [x, y] = createXY(N);
    plot(x, y, '*');
    xlabel('x');
    ylabel('y');
    title(['N = ', num2str(N)]);
end


Ns = [10, 20, 40];
Ks = [5, 18, 35];

%%Problem 2
for i = 1 : 3
    N = Ns(i);
    K = Ks(i);
    [x, y, fxLS, fxLSChol] = LSsolve(K, N);
    figure(i+3)
    plot(x1, f(x1), x, fxLS, '-o')
    legend('reference function', 'result of approximation')
    grid on
    xlabel('x')
    ylabel('f(x)')
    title(['The approximation of f(x) for N = ', num2str(N), ' and K = ', num2str(K)]);
end

%Problem 3
%Dependence of indicator on N and K
K = [4:40];
N = [6:42];
accuracy = zeros(length(K), length(N));
accuracyInf = zeros(length(K), length(N));

for k = 4 : 40
    for n = k+2 : 42
        [x, y, fxLS, ] = LSsolve(k, n);
        accuracy(k-3, n-k-1) = norm(fxLS - y) / norm(y);
        accuracyInf(k-3, n-k-1) = norm(fxLS - y, Inf) / norm(y, Inf);
    end
end

figure(7)
mesh(K, N, log10(accuracy))
xlabel('K');
ylabel('N');
zlabel('\bf \delta_{2}');
title('The dependence of \delta_{2}(K, N) on N and K');
grid on

figure(8)
mesh(K, N, log10(accuracyInf))
xlabel('K');
ylabel('N');
zlabel('\bf \delta_{\infty}');
title('The dependence of \delta_{\infty}(K, N) on N and K');
grid on


%Problem 4
err = [0.01, 0.001, 0.0001];
accuracy = zeros(length(K), length(N));
accuracyInf = zeros(length(K), length(N));
variance = zeros(length(K), length(N));
for i = 1 : 3
    lvDisturbance = err(i);
    for k = 4 : 40
        for n = k+2 : 42
            [x, y, fxLS, ] = LSsolveErr(k, n, lvDisturbance);
            accuracy(k-3, n-k-1) = norm(fxLS - y) / norm(y);
            accuracyInf(k-3, n-k-1) = norm(fxLS - y, Inf) / norm(y, Inf);
            variance(k-3, n-k-1) = var(y);
        end
    end

    figure(2*i + 7)
    mesh(K, N, log10(accuracy))
    xlabel('K');
    ylabel('N');
    zlabel('\bf \delta_{2}');
    title(['The dependence of \delta_{2}(K, N) on \sigma_{y}^{2} with level of disturbance is ', num2str(lvDisturbance)]);
    grid on

    figure(2*i + 8)
    mesh(K, N, log10(accuracyInf))
    xlabel('K');
    ylabel('N');
    zlabel('\bf \delta_{\infty}');
    title(['The dependence of \delta_{\infty}(K, N) on \sigma_{y}^{2}with level of disturbance is ', num2str(lvDisturbance)]);
    grid on
end

%%Function to solve approximation
function[x, y, fxLS, fxLSChol] = LSsolve(K, N)
    [x, y] = createXY(N);
    y = y';
    P = createBase(K, N, x);
    res = (P' * P) \ (P' * y);
    resCB = solveCB(P' * P, P' * y);
    fxLS = P * res;
    fxLSChol = P * resCB;
end


%%Function to solve approximation with error
function[x, y, fxLS, fxLSChol] = LSsolveErr(K, N, err)
    [x, y] = createXY(N);
    y = y';
    yErr = randn(N,1) * err;
    y = y .* (1 + yErr);
    P = createBase(K, N, x);
    res = (P' * P) \ (P' * y);
    resCB = solveCB(P' * P, P' * y);
    fxLS = P * res;
    fxLSChol = P * resCB;
end


%%Function to create data {(xn, yn)|n = 1, ..., N}-------------
function [X, Y] = createXY(N)
    f = @(x) -sin(pi*x).*exp(-x);
    for n = 1 : N
        X(n) = -1 + 2*(n-1)/(N-1);
    end
    Y = f(X);
end


%%Function to creaatebase---------------------------------------
function P = createBase(K, N, x)
    for n = 1 : N
        P(n, 1) = 1;
        P(n, 2) = x(n);
        for j = 3 : K+1
            k = j-1;
            P(n, j) = (2*k-1)/k * x(n) * P(n, j-1) - (k-1)/k * P(n, j-2);
        end
    end
end


%%Cholesky function----------------------------------------------    
function [L] = Cholesky(A)
    N = length(A);
    L = A-A;
    for i = 1  : N
        L(i, i) = sqrt( A(i, i) - L(i, :)*L(i, :)' );

        for j = (i + 1) : N
            L(j, i) = ( A(j, i) - L(i, :)*L(j, :)' )/L(i, i);
        end
    end
end


%%Function to sovle linear system---------------------------------------
function [X] = solveCB(A, b)
    L = Cholesky(A);
    Lt = L';
    [n , ~] = size(A);
    y = zeros(n, 1);
    X = zeros(n, 1);
    
    y(1) = b(1)/L(1, 1);
    for i = 2 : n
        y(i) = (b(i) - L(i, :)*y)/L(i, i);
    end
    
    X(n) = y(n)/Lt(n, n);
    for i = n-1 : -1 : 1
        X(i) = (y(i) - Lt(i, :)*X)/L(i, i);
    end
end