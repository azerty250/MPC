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

N=5;
iter=50;
Q=eye(2);
R=1; %Weight
P = dlyap(A,Q);
r=1;  %Comand

X=zeros(2,iter);
d=0.2;
u=zeros(1,iter);
y=zeros(1,iter);

Xhat=zeros(2,iter+1);
%uhat=zeros(1,iter);
dhat=zeros(1,iter+1);


X(:,1)=[1;2];
y(1)=C*X(:,1)+d;
Xhat(:,1)=[3;0];
dhat(1)=0;

for i=1:iter  %Slide 6-56
    %Obtain Xs and Us
    Xss=sdpvar(2,1);
    Uss=sdpvar(1,1);
    rVar=sdpvar(1,1);
    dhatVar=sdpvar(1,1);

    obj = Uss*Uss;
    con=[];
    con=[con, (eye(2)-A)*Xss==B*Uss];  %Steady state  Bd=0
    con=[con, C*Xss==r - dhatVar];  %Matching the comand
    con=[con, [1;-1]*Uss <= [3;3]];

    opti = optimizer(con, obj, [], dhatVar, [Xss;Uss]);  %ctrl(r) give [x_1; x_2; u] at steady state
    
    vector=opti(dhat(i));
    Xss=vector(1:2);
    Uss=vector(3);
    %Solve the MPC problem
    Xvar=sdpvar(2,N);
    Uvar=sdpvar(1,N);
    obj=0;
    con=[];
    for j=1:N-1
        con=[con, Xvar(:,j+1)==A*Xvar(:,j)+B*Uvar(j)];  %the system is respected
        con=[con,[1;-1]*Uvar(j)<= [3;3] ];          %Constraint on the input
        obj=obj+(Xvar(:,j)-Xss)'*Q*(Xvar(:,j)-Xss)+(Uvar(j)-Uss)'*R*(Uvar(j)-Uss);
    end
    obj=obj+(Xvar(N)-Xss)'*P*(Xvar(N)-Xss);
    ctrl= optimizer(con, obj, [],Xvar(:,1) , Uvar(1));
    
    u(i)=ctrl(Xhat(:,i));
    
    %We have a controller, the state evoluate:
    X(:,i+1)=A*X(:,i)+B*u(i);
    y(i+1)=C*X(:,i+1)+d;
    
    %Estimation of the stade and distrubance
    vector=[A Bd; zeros(1,2),1]*[Xhat(:,i);dhat(i)]+[B;0]*u(i)+L*(C*Xhat(:,i)+Cd*dhat(i)-y(i));
    Xhat(:,i+1)=vector(1:2);
    dhat(i+1)=vector(3);
    disp(i) %To make us keep waiting ;)
end