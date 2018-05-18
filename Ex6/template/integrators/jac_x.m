function [A] = jac_x(X0,U0,f)
% Write your jacobian function here using finite differences

d = [1;1]*0.00001;

A = zeros(length(X0));

A(:,1) = 0.5*(f(X0+d(1),U0)-f(X0-d(1),U0))/d(1);
A(:,2) = 0.5*(f(X0+d(2),U0)-f(X0-d(2),U0))/d(2);


end

