close all;
clear all;
yalmip('clear')

clc;

%%
A = [0.7115 -0.4345; 0.4345 0.8853];
B = [0.2173; 0.0573];
C = [0,1];
Cd = 1;
Bd = zeros(2,1);

N =5;
Q = eye(2);
R = 1;

x0_hat = [3;0];
x0 = [1;2];
d = 0.2;
d0_hat = 0;
u = 0;

%Constraints on u
G = [1; -1];
g = [3;3];

% % Define optimization variables
% x = sdpvar(2,N,'full');
% u = sdpvar(1,N,'full');

A_hat = [A Bd; zeros(1,2) 1]';
B_hat = [C Cd]';
F = [0.5,0.6,0.7];


P = dlyap(A,Q);
L = -place(A_hat,B_hat,F)';

x_hat = x0_hat;
d_hat = d0_hat;
out_hat = [x_hat;d_hat];
x = x0;

for i = 1:100
    x(:,i+1) = A*x(:,i) + B*u;
    y(i) = C*x(:,i)+Cd*d;
    out_hat(:,i+1) = A_hat'*out_hat(:,i) + [B;0]*u + L*(C*x_hat + Cd*d_hat - y(i)); 
    x_hat = out_hat(1:2,i+1);
    d_hat = out_hat(end,i+1);
end

%%
sys = LTISystem('A',A,'B',B);

sys.u.min = -3;
sys.u.max =  3;

sys.x.penalty = QuadFunction(Q);
sys.u.penalty = QuadFunction(R);

Qf = sys.LQRPenalty.weight();
Xf = sys.LQRSet();

Ff = Xf.A;
ff = Xf.b;

%%
% Define constraints and objective
D = eye(2)-A;
M = [D,-B;C,0];

steady_state = M\[zeros(2,1);r];

xs = steady_state(1:2);
us = steady_state(end);

dx = x - xs;


con = [];
obj = 0;
for i = 1:N-1
    con = con + (G*u(:,i) <= g); % Input constraints
    obj = obj + x(:,i)'*Q*x(:,i) + u(:,i)'*R*u(:,i); % Cost function
end

con = con + (Ff*x(:,N) <= ff); % Terminal constraint
obj = obj + dx'*P*dx; % Terminal weight

ops = sdpsettings('solver','quadprog');
% Compile the matrices
ctrl = optimizer(con, obj,ops, x(:,1), u(:,1));
% Can now compute the optimal control input using
[uopt(1), isfeasible] = ctrl{x0};
% isfeasible == 1 if the problem was solved successfully

x_plot = x0;
i =1;

while i < 30
    
    [uopt(i), isfeasible] = ctrl{x_plot(:,i)};

    x_plot(:,i+1) = A*x_plot(:,i) + B*uopt(i); 
      
    i = i + 1;
end

