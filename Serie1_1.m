%% Ex1

A = [4/3 -2/3 ; 1 0];
B = [1 ;0];
C = [-2/3 1];

Q = C'*C + 0.001*eye(2);

R = 0.001;

%% Ex1.2

N = 7;
x0 = [10;10];

H_old = Q;
K = [];

for i=N:-1:1
    
    K_new = -inv(R + B'*H_old*B)*B'*H_old*A;
    K = [K K_new'];
    K_old = K_new;
    H_new = Q + K_old'*R*K_old + (A+B*K_old)'*H_old*(A+B*K_old);
    H_old = H_new;
    
end

K0 = K_new;
x_old = x0;

x = [x0];

for i=N:-1:1
    
    x_new = A*x_old + B*K0*x_old;
    x = [x x_new];
    x_old = x_new;
end

%% Prediction plot
x_pred = zeros(2,N);
x_pred(:,1) = x0;

for i = 1:length(K)
   x_pred(:,i+1) = A*x_pred(:,i) + B.*K(:,i).*x_pred(:,i);
end

plot(x_pred(1,:), x_pred(2,:),'-x');
title('Prediction')
    
%% Ex3

[K_inf,S,e] = dlqr(A,B,Q,R);
cost = 0;
x = x0;
K_cost = -K_inf;

for i = 1:1000

    cost = cost + x'*Q*x + x'*K_cost'*R*K_cost*x;
    
    x = A*x + B*K_cost*x;
    
end

