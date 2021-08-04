clear;
format long; format compact

max_iter = 10;
tol = 1e-6;

% test diagonal matrix
%A = diag([-1,2,1,7,3])
%x = rand(size(A,1),1)

A = [-2 1 4; 1 1 1; 4 1 -2];
x = [1, 2, 1]';

% X is the matrix of eigenvector guesses
[X, eigenvalues] = PowerMethod(A,x,max_iter,tol);

% print the kth eigenvalue/eigenvector iteration
k = 6;
computed_eigenvector = X(:,k)
computed_eigenvalue = eigenvalues(max_iter)


[B,D] = eig(A);
true_eigenvectors = B
true_eigenvalues = D



function [X, rq] = PowerMethod(A,x,iter,tol)
    n = length(x);
    X = zeros(n, iter);  % stores nx1 eigenvectors for each iteration
    rq = zeros(iter, 1); % stores the guess of eigenvalue at each iteration
    q = x;
    q = q/norm(q,2);
    
    for k = 1:iter
        X(:,k) = q;     % save this guess of the eigenvector
        q = A*q;        % compute q_k+1 and override value of q_k
        q = q/norm(q,2); % normalize q_k+1
        
        rq(k) = (q')*A*q;  % compute raleigh quotient
        
        % compute residual and terminate if within tolerance
        if(norm(q-rq(k)*q,2) < tol)
            % remove 0 columns in X if we terminate early
            X(:,all(X == 0)) = [];  % this line is taken from matlab forum
            return;
        end
    end
end