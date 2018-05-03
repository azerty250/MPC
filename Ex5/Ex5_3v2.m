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

A_hat = [A Bd; zeros(1,2) 1]';
B_hat = [C Cd]';
F = [0.6,0.7,0.8];

L = -place(A_hat,B_hat,F)';

x_hat = x0_hat;
d_hat = d0_hat;
out_hat = [x_hat;d_hat];

x_true = x0;
u_est = 0;

D = eye(2)-A;
M = [D,-B;C,0];
temp =  r-Cd*d_hat;

%y(1)=C*x_true(:,1)+d0;

for i = 1:50
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Compute the steady state

    us = sdpvar(1,1,'full');
    xs = sdpvar(2,1,'full');
    r = sdpvar(1);
    d = sdpvar(1);

    % Define constraints and objective for steady state
    con_ss = [];
    obj_ss = 0;

    con_ss = con_ss + (xs == (A*xs+B*us));
    con_ss = con_ss + (C*xs + d == r);
    con_ss = con_ss + (-3<=us<=3);

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
    x = sdpvar(2,N);
    u = sdpvar(1,N);
    
    % Define constraints and objective for mpc
    con_mpc = [];
    obj_mpc = 0;

    for j = 1:N-1

        con_mpc = con_mpc + (x(:,j+1)== A*x(:,j) + B*u(j));
        con_mpc = con_mpc +  (-3<=u<=3);
        obj_mpc = obj_mpc + (x(:,j)-xs)'*Q*(x(:,j)-xs) + (u(j)-us)'*R*(u(j)-us);

    end
    V = (x(:,N)-xs)'*P*(x(:,N)-xs);
    obj_mpc = obj_mpc + V;
    
    ops = sdpsettings('solver','quadprog');
    % Compile the matrices
    ctrl_mpc = optimizer(con_mpc, obj_mpc,ops,x(:,1),u(1));

    u_est(i) = ctrl_mpc(x_true(:,i));
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Compute the estimates
    x_true(:,i+1) = A*x_true(:,i) + B*u_est(i) ;
    y(i) = C*x_true(:,i)+Cd*d0;
    
    out_hat(:,i+1) = A_hat'*out_hat(:,i) + [B;0]*u_est(i) + L*(C*x_hat + Cd*d_hat - y(i)); 
    x_hat = out_hat(1:2,i+1);
    d_hat = out_hat(end,i+1);
end
 %%

figure(1)
plot(out_hat(1,:),out_hat(2,:)),hold on
plot(x_true(1,:),x_true(2,:)),hold off
grid on

legend('Estimated state','True state');
title('Estimate converging to the true value');
xlabel('Position');
ylabel('Speed');

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'converge_final','-dpdf','-r0')

figure(2)
plot(y)
grid on

legend('Achieved output');
title('Output converges to reference');
xlabel('Index');
ylabel('Output y');

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'output','-dpdf','-r0')

figure(3)
plot(u_est)
hold
plot(0:size(u_est,2)-1, 3*ones(size(u_est,2),1), 'r');
plot(0:size(u_est,2)-1, -3*ones(size(u_est,2),1), 'r');
grid on

ylim([-3.25 3.25]);

legend('Computed input');
title('Input does not violate constraints');
xlabel('Index');
ylabel('Input u');

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'input','-dpdf','-r0')
