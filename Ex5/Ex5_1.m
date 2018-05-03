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
N = 100;
for i = 1:N
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

legend('Estimated state','True state');
title('Estimate converging to the true value');
xlabel('State first component');
ylabel('State second component');

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'converge','-dpdf','-r0')

figure(3)
plot(out_hat(3,:))


