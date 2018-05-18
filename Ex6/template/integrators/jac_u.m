function [B] = jac_u(X0,U0,f)
% Write your jacobian function here using finite differences

d = [1;1]*0.00001;

B = zeros(length(X0));

B(:,1) = 0.5*(f(X0,U0+d(1))-f(X0,U0-d(1)))/d(1);
B(:,2) = 0.5*(f(X0,U0+d(2))-f(X0,U0-d(2)))/d(2);


end