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


x0_hat = [3;0];
x0 = [1;2];
d0 = 0.2;
d0_hat = 0;
u = 0;

N = 5;
Q = eye(2);
R = 1;

u = 0;

r = 1;

%Constraints on u
G = [1; -1];
g = [3;3];


A_hat = [A Bd; zeros(1,2) 1]';
B_hat = [C Cd]';
F = [0.6,0.7,0.8];

L = -place(A_hat,B_hat,F)';

x_hat = x0_hat;
d_hat = d0_hat;
out_hat = [x_hat;d_hat];
x = x0;

for i = 1:100
    x(:,i+1) = A*x(:,i) + B*u;
    y(i) = C*x(:,i)+Cd*d0;
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

%% Compute steady state

nu =1;
nx =2;
ny = 1;
D = eye(2)-A;
M = [D,-B;C,0];
temp =  r-Cd*d_hat;

us = sdpvar(1,1,'full');
xs = sdpvar(2,1,'full');
r = sdpvar(1);
d = sdpvar(1);

% Define constraints and objective
con = [];
obj = 0;

% = [con, (G*us(:,i) <= g)]; % Input constraints
%con = [con, (xs(:,i+1) == A*xs(:,i)+B*us(:,i)+Bd*d)];  

con = [con, xs == (A*xs+B*us)];
con = [con, C*xs + d == r];
con = [con, (-3<=us<=3)];

obj = us.*us; % Cost function

ops = sdpsettings('solver','quadprog');

input = [r,d];
output = [xs;us];

% Compile the matrices
ctrl = optimizer(con, obj,ops, input, output);


d = d_hat;
r = 1;
input = [r,d];

% Can now compute the optimal control input using
steady = ctrl{input};


%% 
P = dlyap(A,Q);

x_var = sdpvar(2,N,'full');
u_var = sdpvar(1,N,'full');


V = 0;
for i = 1:N-1
    dx(:,end+1) = A*dx(:,end)+B*du;
    V = V + dx(:,end+1)'*Q*dx(:,end+1) + du'*R*du;
end

% V = V + dx()
