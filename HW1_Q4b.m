clear
format long; format compact

% Intersection of ellipse and line

f = @(x) 4.81*x^2-10.8*x;  % the function, in x instead of t

x_exact=(10.8/4.81)     % the exact value from analytic calculation
intersection=[x_exact,2-0.3*x_exact]  % the exact intersection
tol=1e-6;
n_maxiter=8;
x_secant=zeros(n_maxiter,1);

x_old=1; % We know t=0 is a trivial solution to this problem, start at x0=1
f_old=f(x_old);
x=4; % The largest x-value a point on the ellipse has is 3, so 4 is reasonable for x1
n_iter=n_maxiter;
for k=1:n_maxiter
   f_x = f(x);
   
   if(abs(f_x)<tol) % termination criterion
      n_iter=k-1
      break
   end    
   
   x_new = x-f_x*(x-x_old)/(f_x-f_old);
   
   f_old=f_x;
   x_old=x;   
   x=x_new;
   
   x_secant(k)=x;

end   
x_secant(1:n_iter)
error = x_exact - x