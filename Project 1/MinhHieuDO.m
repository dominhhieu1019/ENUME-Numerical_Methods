clear all
close all
f = @(x) 2/x - 2;

Ns = [3, 10, 20];

for n = 1 : 3
    N = Ns(n);
    X = [1:N]';
    %%Problem 1-----------------------------------------------

    an = 1;
    %problem 2----------------------------------------------
    d = 101;
    a = (an-0.01 : 0.02/d : an+0.01);
    
    figure (1)
    for i = 1 : d+1 
        A = createMatrix(N, a(i));
        detA(i) = det(A);
        
        condA(i) = cond(A);
        
        B = A * X;
        Xn = solveCB(A, B);
        XnMatlab = A\B;
        
        accuracy(i) = norm(Xn - X)/norm(Xn);
        accuracyInf(i) = norm(Xn - X, Inf)/norm(Xn, Inf);
        
        accuracyMatlab(i) = norm(XnMatlab - X)/norm(XnMatlab);
        accuracyInfMatlab(i) = norm(XnMatlab - X, Inf)/norm(XnMatlab, Inf);   
    end
    
    
    figure (1)
    semilogy(a, detA);
    title("The dependence of the determinant on a")
    xlabel("a")
    ylabel("Determinant")
    hold on
    grid on


    figure (2)
    semilogy(a, condA);
    hold on
    grid on
    title('The dependence of the condition number on a')
    xlabel('a')
    ylabel('Condition number')
    
    
    %Problem 4, 5----------------------------------------
    
%     AA = subs(A, x, an);
    figNum = 2*n;
    figure (figNum + 1)
    semilogy(a, accuracy, a, accuracyMatlab);
    legend('CB', 'Matlab')
    hold on
    grid on
    title(['The dependence of accuracy indicator 2 on a (N = ', num2str(N), ')'])
    xlabel('a')
    ylabel('Accuracy indicator 2')
    
    figure (figNum + 2)
    semilogy(a, accuracyInf, a, accuracyInfMatlab);
    legend('CB', 'Matlab')
    hold on
    grid on
    title(['The dependence of accuracy indicator infinitive on a(N = ', num2str(N), ')'])
    xlabel('a')
    ylabel('Accuracy indicator infinitive')
    
%     figure (5)
%     semilogy(a, accuracyMatlab);
%     hold on
%     grid on
%     title("The dependence of accuracy indicator 2 on a (Matlab)")
%     xlabel("a")
%     ylabel("Accuracy indicator 2")
%     
%     figure (6)
%     semilogy(a, accuracyInfMatlab);
%     hold on
%     grid on
%     title("The dependence of accuracy indicator infinitive on a (Matlab)")
%     xlabel("a")
%     ylabel("Accuracy indicator infinitive")
    
end



%%Function to create matrix A ----------------------------------
function A = createMatrix(N, x)
    f = @(x) 2/x - 2;
    for i = 1 : N
        for j = 1 : N
            if i == 1 && j == 1
                A(i,j) = f(x).^2;
            elseif i == 1 || j == 1
                A(i,j) = (-1).^abs(i-j).*2.*f(x);
            else
                A(i,j) = (-1).^(i+j).*4.*min(i,j);
            end
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

%%-----------------------



