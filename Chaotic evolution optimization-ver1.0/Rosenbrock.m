function f = Rosenbrock(x)
% if size(x,1) < 2 error('dimension must be greater one'); end
% -30 <= x <= 30 
  f = 100*sum((x(1:end-1).^2 - x(2:end)).^2) + sum((x(1:end-1)-1).^2);