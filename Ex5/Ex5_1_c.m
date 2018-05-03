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
F = [0.6,0.7,0.8];

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

figure(1)
plot(out_hat(1,:),out_hat(2,:)),hold on
plot(x(1,:),x(2,:)),hold off
grid on

figure(3)
plot(out_hat(3,:))


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
con = [];
obj = 0;
for i = 1:N-1
    con = con + (G*u(:,i) <= g); % Input constraints
    obj = obj + x(:,i)'*Q*x(:,i) + u(:,i)'*R*u(:,i); % Cost function
end

con = con + (Ff*x(:,N) <= ff); % Terminal constraint
obj = obj + x(:,N)'*Qf*x(:,N); % Terminal weight

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

%%
figure
plot(0:size(x_plot,2)-1,x_plot(1,:));
hold
plot(0:size(x_plot,2)-1, 5*ones(size(x_plot,2),1), 'r');
plot(0:size(x_plot,2)-1, -5*ones(size(x_plot,2),1), 'r');
xlim([0 25]);
ylim([-5.25 5.25]);
legend('Position','Constraints');
title('Plot of the position');
xlabel('Time (step)');
ylabel('Computed value');

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'4_2_pos','-dpdf','-r0')

figure   
plot(0:size(x_plot,2)-1,x_plot(2,:));
hold
plot(0:size(x_plot,2)-1, 0.2*ones(size(x_plot,2),1), 'r');
plot(0:size(x_plot,2)-1, -0.2*ones(size(x_plot,2),1), 'r');
xlim([0 25]);
ylim([-0.25 0.25]);
legend('Speed','Constraints');
title('Plot of the speed');
xlabel('Time (step)');
ylabel('Computed value');

clear('set');
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'4_2_vel','-dpdf','-r0')

figure
plot(0:size(uopt,2)-1,uopt(:));
hold
plot(0:size(x_plot,2)-1, 1.75*ones(size(x_plot,2),1), 'r');
plot(0:size(x_plot,2)-1, -1.75*ones(size(x_plot,2),1), 'r');
xlim([0 25]);
ylim([-2 2]);
legend('Position','Constraints');
title('Plot of the input');
xlabel('Time (step)');
ylabel('Computed value');

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
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
F = [0.95,0.96,0.97];

L = place(A_hat,B_hat,F)';

x_hat = x0_hat;
d_hat = d0_hat;
out_hat = [x_hat;d_hat];

for i = 1:100
    y(i) = C*x_hat+Cd*d;
    out_hat(:,i+1) = A_hat'*out_hat(:,i) + [B;0]*u + L*(C*x_hat + Cd*d_hat - y(i)) 
    x_hat = out_hat(1:2,i+1);
    d_hat = out_hat(end,i+1);
end

figure(1)
plot(out_hat(1,:),out_hat(2,:))
grid on

figure(2)
plot(out_hat(3,:))



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
con = [];
obj = 0;
for i = 1:N-1
    con = con + (G*u(:,i) <= g); % Input constraints
    obj = obj + x(:,i)'*Q*x(:,i) + u(:,i)'*R*u(:,i); % Cost function
end

con = con + (Ff*x(:,N) <= ff); % Terminal constraint
obj = obj + x(:,N)'*Qf*x(:,N); % Terminal weight

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

%%
figure
plot(0:size(x_plot,2)-1,x_plot(1,:));
hold
plot(0:size(x_plot,2)-1, 5*ones(size(x_plot,2),1), 'r');
plot(0:size(x_plot,2)-1, -5*ones(size(x_plot,2),1), 'r');
xlim([0 25]);
ylim([-5.25 5.25]);
legend('Position','Constraints');
title('Plot of the position');
xlabel('Time (step)');
ylabel('Computed value');

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'4_2_pos','-dpdf','-r0')

figure   
plot(0:size(x_plot,2)-1,x_plot(2,:));
hold
plot(0:size(x_plot,2)-1, 0.2*ones(size(x_plot,2),1), 'r');
plot(0:size(x_plot,2)-1, -0.2*ones(size(x_plot,2),1), 'r');
xlim([0 25]);
ylim([-0.25 0.25]);
legend('Speed','Constraints');
title('Plot of the speed');
xlabel('Time (step)');
ylabel('Computed value');

clear('set');
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'4_2_vel','-dpdf','-r0')

figure
plot(0:size(uopt,2)-1,uopt(:));
hold
plot(0:size(x_plot,2)-1, 1.75*ones(size(x_plot,2),1), 'r');
plot(0:size(x_plot,2)-1, -1.75*ones(size(x_plot,2),1), 'r');
xlim([0 25]);
ylim([-2 2]);
legend('Position','Constraints');
title('Plot of the input');
xlabel('Time (step)');
ylabel('Computed value');

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'4_2_input','-dpdf','-r0')