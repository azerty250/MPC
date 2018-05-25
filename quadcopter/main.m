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

N = 20; 
T = 2;

angleMax = deg2rad(10);
vAngleMax = deg2rad(15);

%Use for stage cost --> to tune
Q = 10*eye(7); 
R = eye(4);

A1 = sys.A;
B1 = sys.B;

% Terminal controller
[K1, P1] = dlqr(A1,B1,Q,R); K1 = -K1;

x0 = [-1 deg2rad(10) deg2rad(-10) deg2rad(120) 0 0 0];

x = sdpvar(7, N);
u = sdpvar(4,N-1);
constraints = [];
objective = 0;

for i = 1:N-1
    % Dynamics
    constraints = [constraints, x(:,i+1) == A1*x(:,i) + B1*u(:,i)];
    % Input constraints
    constraints = constraints + (0 <= u(:,i) <= 1);    
    % Cost
    objective = objective + u(:,i)'*R*u(:,i);
end

for i = 2:N
    %constraints on alpha, beta, alpha dot and beta dot
    constraints = constraints + (-angleMax <= x(2:3,i) <= angleMax);
    constraints = constraints + (-vAngleMax <= x(5:6,i) <= vAngleMax);
    
    %Objective
    objective = objective + (x(:,i))'*Q*(x(:,i));

end

objective = objective + (x(:,N))'*P1*(x(:,N));

options = sdpsettings('solver','sedumi');
innerController = optimizer(constraints, objective, options, x(:,1), u(:,1));
simQuad( sys, innerController, x0, T);

pause

%% Reference tracking - no disturbance, no invariant sets
fprintf('PART II - reference tracking...\n')

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






