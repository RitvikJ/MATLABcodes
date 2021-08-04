clc;
clear;
format long; format compact;

% f = 1 for x>= 0 and -1 for x < 0
f = @(x) (heaviside(x)-heaviside(-x))+(1-abs(sign(x)));

% lines 9-15,29-32 based on code written by kiara chen in my worksheet 5 group
n = 8;
a = -1;
b = 1;
x = linspace(a,b,20*(b-a));
num_points = 2.^[3:8];

[nodes,p] = InterpolatingPolynomial(f,n,a,b);


% plot(nodes,f(nodes),'ko'); hold on;
% L1 = plot(x,f(x),'k-');
% L2 = plot(x,polyval(p,x),'r--');
% legend([L1,L2],'f(x)','p(x)');
% title('Degree =',n);
% hold off;

PlotError(f,x,a,b,num_points)

function [nodes,p_n] = InterpolatingPolynomial(f,n,a,b)
    
    k = [0:n];  % code from worksheet, creates chebyshev nodes
    nodes = cos((2*k+1) * pi / (2 * (n+1)));
    nodes = (nodes + (a + b)/2) * ((b - a)/2);
    
    % create a matrix to store the coefficients of each L_k
    poly_coeff = zeros(n+1,n+1)';

    for i = 1:n+1
        prod = 1;
        for j = setdiff(1:n+1,i)  % all 1,2..n except i
            prod = prod*(nodes(i)-nodes(j)); % denominator of L_i
        end
        % stores coefficients of L_i in the ith column
        poly_coeff(i,:) = f(nodes(i))*(1/prod)*(poly(nodes(setdiff(1:n+1,i))));
    end

    % interpolating polynomial is the sum of lagrange polynomials
    p_n = sum(poly_coeff);  
end


function PlotError(f,x,a,b,points)
    for j = 1:length(points)
        % interpolate for all n given in the vector 'points'
        [~,p] = InterpolatingPolynomial(f,points(j),a,b);
        semilogy(x,abs(f(x)-polyval(p,x))); hold on;
    end
    title('Log Error Plot');
end