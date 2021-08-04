clear;
format long; format compact

U = [1 2 6 -1; 0 3 1 0; 0 0 4 -1; 0 0 0 2];  % initialize matrix
b = [-1; -3; -2; -4];  % initialize vector b
x_sol = backward(U,b)  % print solution
Ux_minus_b = checkSol(U,x_sol,b) % if Ux-b=0, the solution is valid 

function x = backward(A,b)
    % validate input
    len = length(b);
    [n1 n2] = size(A);
    if(n1 ~= n2)
        error('Error, the matrix is not square')
    elseif(n1 ~= len)
        error('Error, the matrix and vector sizes do not match')
    end
    
    x = zeros(len,1);  % empty vector of solutions x
    
    for k=len:-1:1
        % solve upwards from x_n to x_1
        % the value of x_n is simply b_n divided by the nth pivot
        % and stored in x
        x(k) = b(k)/A(k,k);
        
        % update vector b in place as each x_k is found
        b = b - (A(1:len,k).*x(k));        
    end
    
end

function y = checkSol(U,x,b)
    y = U*x - b;
end