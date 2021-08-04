clear;
format long; format compact

n = 6  % number of points, including endpoints
a = 0;   % left endpoint
b = 1;   % right endpoint

f = @(x) (exp(x).*sin(x)).^2;
f_value = integral(f,a,b);

area = Simpson(f,n,a,b)
error = abs(f_value-area)/f_value

gauss_x = [-sqrt(3/5),0,sqrt(3/5)];
gauss_w = [5/9,8/9,5/9];
gauss_estimate = sum(f(gauss_x).*gauss_w)
gauss_error = abs(f_value - gauss_estimate)/f_value



function area = Simpson(f,n,a,b)
    h = (b-a)/(n-1);  % interval length
    nodes = [a:h:b];         % nodes x_k
    
    % this has very fast convergence
    midnodes = zeros(n-1,1); % extra nodes between each given node, x_k+1/2
    for k = 1:length(midnodes)
        midnodes(k) = 0.5*(nodes(k)+nodes(k+1));
    end
    
    sum1 = (f(nodes(1))+f(nodes(n)))/6;
    sum2 = sum(f(nodes(2:n-1)))/3;
    sum3 = sum(f(midnodes))*(2/3);
    
    area = h*(sum1+sum2+sum3);
    
    %-------------------------
    % first attempt
    % MUCH MUCH  slower convergence with this code
%     sum = (nodes(1)+nodes(n))/3;
%     
%     for k = 2:length(nodes)-1
%         if mod(k,2)==0
%             sum = sum + (4/3)*f(nodes(k));
%         else
%             sum = sum + (2/3)*f(nodes(k));
%         end
%     end
%     
%     area = h*sum;
end