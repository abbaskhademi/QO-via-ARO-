function [LB,UB,x,Time]=MILO_QO(Q,c,A,b)
% Code written by Abbas Khademi 
% Requirements: Yalmip & Gurobi

%______QO Problem Definitio______%
% min    x'Qx + c'x              %
% s.t.   Ax = b, x >= 0.         %
%________________________________%

%=================================================================%
% ____________ INPUT:                                             %
% Q ........... Matrix of the quadratic term in the objective     %
% c ........... Vector of the linear term in the objective        %
% A ........... Coefficient matrix for equality constraints       %
% b ........... Right-hand side vector for equality constraints   %
% ____________ OUTPUT:                                            %
% UB ........... Upper bound value                                %
% LB ........... Lower bound value                                %
% x ............  solution vector                                 %
% Time ......... Total solver time                                %
%=================================================================%
% Dimension
[m,n] = size(A);  

%% bounds on primal variables
x = sdpvar(n, 1); 
U=[];                   
Time_big_M=0;
for i=1:n
    sol = optimize([A*x == b, x >= 0], -x(i), sdpsettings('verbose', 0, 'solver', 'gurobi'));
    Time_big_M = Time_big_M+sol.solvertime;
    U = [U; value(x(i))];
end
%% bounds on the dual variables
yalmip clear;
% Decision variables
x = sdpvar(n,1);
lambda = sdpvar(n,1);
mu = sdpvar(length(b),1);
rho = sdpvar(n,1);
X = sdpvar(n,n,'symmetric');
V=[];
TT=U*U';
for i=1:n
    % Objective function
    Objective = -lambda(i); 

    % Constraints
    Constraints = [2*Q*x + c + A'*mu - lambda + rho == 0, ...
               trace(2*Q*X) + c'*x + b'*mu + U'*rho == 0, ...
                0 <= X(:) <= kron(U, U)', ...
                0 <= x <= U, ...
                lambda >= 0, ...
                rho >= 0];
    % Solve the problem
    options = sdpsettings('solver', 'gurobi','verbose', 0); 
    sol = optimize(Constraints, Objective, options);
    V = [V; value(lambda(i))];
    Time_big_M=Time_big_M+sol.solvertime;
end

%%    MILO-QO
% Define the variables
yalmip clear;
x = sdpvar(n, 1); 
lambda = sdpvar(n, 1); 
z = binvar(n, 1);
mu = sdpvar(m,1); 
% Constraints
Constraints = [];
Constraints = [Constraints, A*x == b, x >= 0, lambda >= 0, 2*Q*x+c+A'*mu-lambda==0];
Constraints = [Constraints, lambda<=(1-z).*V, x<=z.*U];
% Objective
Objective = 0.5*(c'*x-b'*mu);
options = sdpsettings('verbose', 0, 'solver', 'gurobi','gurobi.TimeLimit',3000-Time_big_M,'savesolveroutput',1);
sol9 = optimize(Constraints, Objective, options);
disp(yalmiperror(sol9.problem))
disp('____________MILO-QO__________'); 
%fprintf('Time Parameter: %s\n', mat2str(Time_big_M));
LB= sol9.solveroutput.result.objbound;
fprintf('LB MILP:   %s\n', mat2str(LB));
UB=value(0.5*(c'*x-b'*mu));
fprintf('UB MILP:   %s\n', mat2str(UB));
Time=sol9.solvertime+Time_big_M;
fprintf('Time MILO: %s\n', mat2str(Time));
disp('_____________________________'); 
end