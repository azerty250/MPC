yalmip('clear')
clear all

% Model data
A = [0.9752 1.4544; -0.0327 0.9315];
B = [0.0248; 0.0327];

nx = 2; % Number of states
nu = 1; % Number of inputs

% MPC data
Q = 10*eye(2);
R = 1;
N = 10;
% Initial state
x0 = [3;0];

u = sdpvar(repmat(nu,1,N),repmat(1,1,N));
x = sdpvar(repmat(nx,1,N+1),repmat(1,1,N+1));

constraints = [];
objective = 0;
for k = 1:N
 objective = objective + norm(Q*x{k},1) + norm(R*u{k},1);
 constraints = [constraints, x{k+1} == A*x{k} + B*u{k}];
 constraints = [constraints, -1 <= u{k}<= 1, -5<=x{k+1}<=5];
end

controller = optimizer(constraints, objective,[],x{1},[u{:}]);

x = [3;1];
clf;
hold on
implementedU = [];
for i = 1:15
  U = controller{x};  
  stairs(i:i+length(U)-1,U,'r')
  x = A*x + B*U(1);
  pause(0.05)
  stairs(i:i+length(U)-1,U,'k')
  implementedU = [implementedU;U(1)];
end
stairs(implementedU,'b')