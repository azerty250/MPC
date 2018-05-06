

%define the dynamic functions
gamma = 3;
alpha = 1.6;
beta = 0.01;
U0 = [0, -0.01];
X0 = [0, 0.5];
v_dot = @(v) gamma*U0(1) - alpha*U0(2)*v -beta*(v^3);
lambda_dot = @(v) v;


%Euler forward method
h=1.5;                                             % step size
x = 0:h:3;                                         % Calculates upto y(3)
y = zeros(1,length(x)); 
y(1) = 5;                                          % initial condition
F_xy = @(t,r) 3.*exp(-t)-0.4*r;                    % change the function as you desire

for i=1:(length(x)-1)                              % calculation loop
    k_1 = F_xy(x(i),y(i));
    k_2 = F_xy(x(i)+0.5*h,y(i)+0.5*h*k_1);
    k_3 = F_xy((x(i)+0.5*h),(y(i)+0.5*h*k_2));
    k_4 = F_xy((x(i)+h),(y(i)+k_3*h));

    y(i+1) = y(i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;  % main equation
end