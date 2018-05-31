function [x_next] = RK4_1(X,U,h,f)
%
% Inputs : 
%    X, U current state and input
%    h    sample period
%    f    continuous time dynamics f(x,u)
% Returns
%    State h seconds in the future
%

% Runge-Kutta 4 integration
% write your function here


    k_1 = h*f(X,U);
    k_2 = h*f(X+0.5*k_1,U+0.5*h);
    k_3 = h*f(X+0.5*k_2,U+0.5*h);
    k_4 = h*f(X+k_3,U+h);

    x_next = X + (1/6)*(k_1+2*k_2+2*k_3+k_4);  % main equation
end