
N = 1;

%define the dynamic functions
gamma = 3;
alpha = 1.6;
beta = 0.01;
U0 = [0; -0.01];
X0 = [0; 0.5];

u = U0;
x = X0;

v_dot = @(v,u) gamma*u(1) - alpha*u(2)*v -beta*(v.^3);
lambda_dot = @(v) v;

h=0.1;   % step size
  
f_discrete = X0;

%%
%Euler forward method
for i = 1:N
    
    f_discrete(1,i+1) = f_discrete(1,i) + h*v_dot(f_discrete(2,i), u(:,i)); 
    f_discrete(2,i+1) = f_discrete(1,i) + h*lambda_dot(f_discrete(2,i)); 
    
end

f_discrete

%%

%RK4
for i=1:N                              % calculation loop
    k_1 = v_dot(f_discrete(2,i), u(:,i));
    k_2 = v_dot(f_discrete(2,i)+0.5*h*k_1,u(:,i)+0.5*h);
    k_3 = v_dot((f_discrete(2,i)+0.5*h*k_2),(u(:,i)+0.5*h));
    k_4 = v_dot((f_discrete(2,i)+k_3*h),(u(:,i)+h));

    f_discrete(2,i+1) = f_discrete(2,i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;  % main equation
    
    k_1 = lambda_dot(f_discrete(1,i));
    k_2 = lambda_dot(f_discrete(1,i)+0.5*h*k_1);
    k_3 = lambda_dot((f_discrete(1,i)+0.5*h*k_2));
    k_4 = lambda_dot((f_discrete(1,i)+k_3*h));

    f_discrete(1,i+1) = f_discrete(1,i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;
end

f_discrete