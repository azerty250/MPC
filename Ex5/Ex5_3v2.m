clc;
clear all;

A = [0.7115 -0.4345; 0.4345 0.8853];
B = [0.2173; 0.0573];
C = [0,1];
Cd = 1;
Bd = zeros(2,1);

N = 5;
Q = eye(2);
R = 1;

x0_hat = [3;0];
x0 = [1;2];
d0 = 0.2;
d0_hat = 0;
u = 0;

r = 1;

G = [1; -1];
g = [3;3];

A_hat = [A Bd; zeros(1,2) 1]';
B_hat = [C Cd]';
F = [0.6,0.7,0.8];

L = -place(A_hat,B_hat,F)';

x_hat = x0_hat;
d_hat = d0_hat;
out_hat = [x_hat;d_hat];
x_est = x0;

u_est = 0;

x = sdpvar(2,N,'full');
u = sdpvar(1,N,'full');

for i = 1:50
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Compute the estimates
    x_est(:,i+1) = A*x_est(:,i) + B*u_est(i);
    y(i) = C*x_est(:,i)+Cd*d0;
    out_hat(:,i+1) = A_hat'*out_hat(:,i) + [B;0]*u_est(i) + L*(C*x_hat + Cd*d_hat - y(i)); 
    x_hat = out_hat(1:2,i+1);
    d_hat = out_hat(end,i+1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Compute the steady state
    D = eye(2)-A;
    M = [D,-B;C,0];
    temp =  r-Cd*d_hat;

    us = sdpvar(1,1,'full');
    xs = sdpvar(2,1,'full');
    r = sdpvar(1);
    d = sdpvar(1);

    % Define constraints and objective for steady state
    con_ss = [];
    obj_ss = 0;

    con_ss = [con_ss, xs == (A*xs+B*us)];
    con_ss = [con_ss, C*xs + d == r];
    con_ss = [con_ss, (-3<=us<=3)];

    obj_ss = us.*us; % Cost function

    ops = sdpsettings('solver','quadprog');

    input = [r,d];
    output = [xs;us];

    % Compile the matrices
    ctrl_ss = optimizer(con_ss, obj_ss,ops, input, output);

    d = d_hat;
    r = 1;
    input = [r,d];

    % Can now compute the optimal control input using
    steady = ctrl_ss{input}; 
    
    xs = steady(1:2);
    us = steady(end);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %MPC part 

    P = dlyap(A,Q);

    % Define constraints and objective for mpc
    con_mpc = [];
    obj_mpc = 0;
    
%     sys = LTISystem('A',A,'B',B);
%     sys.u.min = -3;
%     sys.u.max =  3;
% 
%     sys.x.penalty = QuadFunction(Q);
%     sys.u.penalty = QuadFunction(R);
% 
%     Qf = sys.LQRPenalty.weight();
%     Xf = sys.LQRSet();
% 
%     Ff = Xf.A;
%     ff = Xf.b;

    for j = 1:N-1

        con_mpc = [con_mpc, x(:,j+1)== A*x(:,j) + B*u(j)];
        con_mpc = [con_mpc,  (-3<=u<=3)];
        obj_mpc = [obj_mpc, (x(:,j)-xs)'*Q*(x(:,j)-xs) + (u(j)-us)'*R*(u(j)-us)];

    end
    V = (x(:,N)-xs)'*P*(x(:,N)-xs);
    obj_mpc = [obj_mpc,  V];
    
    ops = sdpsettings('solver','quadprog');
    % Compile the matrices
    ctrl_mpc = optimizer(con_mpc, obj_mpc,ops,x,u);
        
end



