% loads:
%    hovering equilibrium (xs,us)
%    continuous time matrices Ac,Bc of the linearization
%    matrices sys.A, sys.B of the inner-loop discretized with sampling period sys.Ts
%    outer controller optimizer instance
load('quadData.mat')
addpath('./cplex');
outerController = getOuterController(Ac, 'sedumi');
disp('Data successfully loaded')



%% %%%%%%%%%%%%%% First MPC controller %%%%%%%%%%%%%%%%%%%

N = 10; 
T = 6;

angleMax = deg2rad(10);
vAngleMax = deg2rad(15);

%State cost --> yaw is less important and can take more time
Q = diag([10,10,10,1,0,0,0]);

%Motor input cost --> all important
R = eye(4);

A1 = sys.A;
B1 = sys.B;

% Terminal controller
[K1, P1] = dlqr(A1,B1,Q,R); 
K1 = -K1;

x0 = [-1 deg2rad(10) deg2rad(-10) deg2rad(120) 0 0 0]';

x1 = sdpvar(7, N); %[zdot alpha beta gamma alphadot betadot gammadot]
u1 = sdpvar(4,N-1); %[mot1 mot2 mot3 mot4]
constraints = [];
objective = 0;

for i = 1:N-1
    % Dynamics
    constraints = constraints + (x1(:,i+1) == A1*x1(:,i) + B1*u1(:,i));
    % Input constraints
    constraints = constraints + (zeros(4,1)-us <= u1(:,i) <= ones(4,1)-us);    
    % Cost
    objective = objective + u1(:,i)'*R*u1(:,i);
end

for i = 2:N
    %constraints on alpha, beta, alpha dot and beta dot
    constraints = constraints + (-angleMax <= x1(2:3,i) <= angleMax);
    constraints = constraints + (-vAngleMax <= x1(5:6,i) <= vAngleMax);
    
    %Objective
    objective = objective + (x1(:,i))'*Q*(x1(:,i));

end

objective = objective + (x1(:,N))'*P1*(x1(:,N));

options = sdpsettings('solver','sedumi');
innerController = optimizer(constraints, objective, options, x1(:,1), u1(:,1));
simQuad(sys,innerController,x0,T);

pause

%% Reference tracking - no disturbance, no invariant sets
fprintf('PART II - reference tracking...\n')

C = [eye(4),zeros(4,3)];

%define a reference of 4 [zdot alpha beta gamma]'
r1 = [0 1 1 1]';

%solve C*xr = r1
xr = C\r1;
ur = B1\(xr - A1*xr);

pause

%% Nonlinear model simulation - no disturbance
fprintf('Running the FIRST NL model simulation...\n')
sim('simulation1.mdl') 

pause

%% Disturbance estimation
%estimator


%% Offset free MPC
fprintf('PART III - OFFSET FREE / Disturbance rejection...\n')

pause

%% Final simulation
fprintf('Running the FINAL NL model simulation...\n')
sim('simulation2.mdl') 
pause
%% BONUS - Slew rate constraints
% run after doing nonlinear simulations otherwise the NL simulations won't
% work (because of the additional controller argument)
fprintf('BONUS - SLEW RATE CONSTRAINTS...\n')






