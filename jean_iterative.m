
% Gauss_seidel and Jacobi

clc; clear; 
format compact;

A = [10 -1 7 -18;2 10 -4 -20;-1.5 6 -20 -2;2 5 -2 30]
b = [2; 1; 9; -18]

% Check for diagonal dominance:
n = size(A);
count = 0;
for k = 1:n
    if abs(A(k,k)) > sum(abs(A(k,:)))-abs(A(k,k))
        count = count+1;
    end
    fprintf("row %d : |diag| %d, |sum| %d\n",k,abs(A(k,k)),sum(abs(A(k,:)))-abs(A(k,k)))
end
if count == n
    disp("The matrix is diagonally dominant")
else
    disp("The matrix is not diagonally dominant")
end


At = A; 
for k = 1:n
    At(k,k) = 0;
end

xold = zeros(n(1),1); xvec = xold; kmax = 100; tol = 10^-7;
count = 2; err = 10; xnew = xold;


while (err>tol) && (count<= kmax)
    for k = 1:n
        xnew(k) = (b(k) - At(k,:)*xnew)/A(k,k); % xold for Jacobi
    end
    
    xvec(:,count) = xnew;
    count = count + 1;
    err = max(abs((xnew-xold)./xnew));
    xold = xnew;
end

d = det(A)
c = cond(A)
r = rank(A)
disp("Build in gaussian")
G = A\b

fprintf("Iterations: %d \n",count-2)

xvec(:,count-1)

if count >= kmax
    disp("Error: method did not converge")
end




        

