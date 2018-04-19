clc
close all
clear all

A = [0.9752 1.4544;
     -0.0327 0.9315];
B = [0.0248;
     0.0327];
H = [1 0;
     -1 0;
     0 1;
     0 -1];
h = [5;5;0.2;0.2];
G = [1;-1];
g = [1.75;1.75];
 
X = Polyhedron(H,h);
U = Polyhedron(G,g);

x0 = [3;0];
%% Q and R initial
Q = 10*eye(2);
R = 10;

% %% Q and R test
% Q = 10*eye(2);
% R = 5;

%%
N = 10;
[X_lqr O_inf] = inv_set_LQR(X,U,A,B,N,Q,R);

plot([X X_lqr O_inf])

figure 
O_inf.plot()
%% MPT3
sys = LTISystem('A',A,'B',B);
sys.x.max = [5;0.2];
sys.x.min = [-5;-0.2];
sys.u.max = 1.75;
sys.u.min = -1.75;


sys.x.penalty = QuadFunction(Q);
sys.u.penalty = QuadFunction(R);

gain = sys.LQRGain;
Qf = sys.LQRPenalty.weight;
myset = sys.LQRSet;

figure
myset.plot()

%% 
F = H; f=h;

H = blkdiag(kron(eye(N-1),Q), Qf, kron(eye(N),R));
h = zeros(size(H,1),1);

Ff = myset.A;
ff = myset.b;

M = G; 
m = g;

G = blkdiag(kron(eye(N-1),F), Ff, kron(eye(N),M));
%g = kron(eye(N-1),f, ff, eye(N), m);
g = [kron(ones(N-1,1),f); ff; kron(ones(N,1),m)];


T =[kron(diag(ones(N-1,1),-1),-A) + kron(eye(N),eye(2)) kron(eye(N),-B)];
t = [A;zeros(2*(N-1),2)];


[zopt, fval, flag] = quadprog(H, h, G, g, T, t*x0);


i = 1;
x_plot = x0;
uopt = [];
while i < 30
    [zopt, fval, flag] = quadprog(H, h, G, g, T, t*x_plot(:,end));
    x_plot(:,end+1) = [zopt(1);zopt(2)];
    uopt(end+1) = [zopt(length(x0)*N+1)];
    i = i+1;
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
print(gcf,'4_1_pos_R10','-dpdf','-r0')

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
print(gcf,'4_1_vel_R10','-dpdf','-r0')

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
print(gcf,'4_1_input_R10','-dpdf','-r0')





