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

N = 5;
Q = eye(2);
R = 1;

u = 0;

r = 1;

%Constraints on u
G = [1; -1];
g = [3;3];


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

%% Compute steady state

D = eye(2)-A;
M = [D,-B;C,0];

steady_state = M\[zeros(2,1);r];

xs = steady_state(1:2);
us = steady_state(end);

% Define constraints and objective
con = [];
obj = 0;

for i = 1:N-1
    con = con + (G*u(:,i) <= g); % Input constraints
    obj = obj + u(:,i).*u(:,i); % Cost function
end

con = con + (Ff*x(:,N) <= ff); % Terminal constraint
obj = obj + x(:,N)'*Qf*x(:,N); % Terminal weight

ops = sdpsettings('solver','quadprog');
% Compile the matrices
ctrl = optimizer(con, obj,ops, x(:,1), u(:,1));
% Can now compute the optimal control input using
[uopt(1), isfeasible] = ctrl{x0};
% isfeasible == 1 if the problem was solved successfully


