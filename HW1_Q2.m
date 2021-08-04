clear
format long; format compact

f = @(x) x^5-3*x^2+1;
g = @(x) 5*x^4-6*x;  %f'(x) for newton's method
tol=1e-1;
tol_newton=1e-12;
n=100;

exact_roots=roots([1 0 0 -3 0 1])
roots_bisec=MultipleBisection(f,-2,2,n,tol)
roots_newton= MultipleNewton(f,g,roots_bisec,n,tol_newton)

function candidates = FindCandidates(f,a,b)
    
    % Find where the function changes sign
    % with max error of 0.25
    candidates=[];
    x_left=f(a);
    for k=-1.75:0.25:2
        x_right=f(k);
        if(x_left*x_right<0)
            candidates(end+1)=x_right;
        end
        x_left = x_right;
    end
    
    candidates(end+1)=b;
    
end

function [x] = SimpleBisection(f,a,b,n,tol)

   f_a=f(a);
   f_b=f(b);
   for k=0:n-1

      x=(a+b)/2; % midpoint of interval

      % Stop as soon as we have reached desired tlerance
      if exist('tol','var') % tol(erance) is an optional argument
         if((b-a)/2<tol)
            n_iter=k; % return how many iterations we took
            break % Already accurate enough
         end
      end      
                     
      f_x=f(x); % Compute f(x) only ONCE per iteration
      
      if(f_x*f_a<0)
         b=x; f_b=f_x;
      else % Should this be elseif(f_x*f_b<0) instead?
         a=x; f_a=f_x;
      end
      
   end
   
   % True answer is somewhere in [a,b]      
   x=(a+b)/2; % Best guess for answer
   
end


function roots = MultipleBisection(f,a,b,n,tol)
    
    % Run bisection on each interval found from FindCandidates()
    roots = zeros(3,1);
    xvals = FindCandidates(f,a,b);
    for k=1:3
        roots(k) = SimpleBisection(f,xvals(k),xvals(k+1),n,tol);
    end
end


%------------------------------------------------------------------
% Functions for Question 2b


function roots = MultipleNewton(f,g,xvals,n,tol_newton)
    
    roots = zeros(length(xvals),1);
    for k=1:length(xvals)
        [r_k,~,~] = Newton(f,g,xvals(k),n,tol_newton);
        roots(k)=r_k; %only stores the value of the root
    end
    
end

function [x,n_iter,x_history] = Newton(f,g,x0,n,tol)
   
   x=x0;  % Initial value
   f_x=f(x);  % f(1.5)
   g_x=g(x);  % f'(1.5)
   n_iter=n;
   for k=1:n
          
      x = x - (f_x/g_x); % Newton with f and g = f'
      
      x_history(k)=x;
            
      f_x=f(x); % Compute f(x) only ONCE per iteration
      g_x=g(x);
      
      if exist('tol','var')
          if(abs(f_x)<tol)
             n_iter=k; % How many iterations we took
             break % Already accurate enough
          end
      end
      
   end
      
end