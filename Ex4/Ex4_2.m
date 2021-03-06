close all;
clear all;
yalmip('clear')

clc;

%%
A = [0.9752 1.4544; -0.0327 0.9315];
B = [0.0248; 0.0327];
N = 10;

x0 = [3;0];

Q = 100*eye(2); %cost we want
R = 1;

H = [1 0; -1 0; 0 1; 0 -1];
h = [5; 5; 0.2; 0.2];

G = [1; -1];
g = [1.75;1.75];

% Define optimization variables
x = sdpvar(2,N,'full');
u = sdpvar(1,N,'full');

sys = LTISystem('A',A,'B',B);

sys.x.min = [ -5 ; -0.2];
sys.x.max = [  5 ;  0.2];
sys.u.min = -1.75;
sys.u.max =  1.75;

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
con = con + (x(:,i+1) == A*x(:,i) + B*u(:,i)); % System dynamics
con = con + (H*x(:,i) <= h); % State constraints
con = con + (G*u(:,i) <= g); % Input constraints
obj = obj + x(:,i)'*Q*x(:,i) + u(:,i)'*R*u(:,i); % Cost function

obj = obj + x(:,i)'*Q*x(:,i) + u(:,i)'*R*u(:,i);
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