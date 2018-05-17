function [x_next] = RK4(X,U,h,f)
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


    k_1 = f(X,U);
    k_2 = f(X+0.5*h*k_1,U+0.5*h);
    k_3 = f(X+0.5*h*k_2),U+0.5*h));
    k_4 = f((X+k_3*h),U+h));

    x_next = X + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;  % main equation
%error('Remove this when you have implemented the function')
