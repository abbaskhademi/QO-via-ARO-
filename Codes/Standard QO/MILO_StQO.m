function [LB,UB,x_value,Time]=MILO_StQO(Q)
% Code written by Abbas Khademi 
% Requirements: Yalmip & Gurobi

%______StQO Problem Definitio____%
% min    x'Qx                    %
% s.t.   e'x = 1, x >= 0.        %
%________________________________%

%=================================================================%
% ____________ INPUT:                                             %
% Q ........... Matrix of the quadratic term in the objective     %
% ____________ OUTPUT:                                            %
% UB ........... Upper bound value                                %
% LB ........... Lower bound value                                %
% x_value ............  solution vector                           %
% Time ......... Total solver time                                %
%=================================================================%
n = size(Q, 1);

% Define the decision variables
x = sdpvar(n, 1);          % Vector x
z = sdpvar(n, 1);          % Vector z
y = binvar(n, 1);          % Binary vector y
alpha = sdpvar(1);         % Scalar alpha
ell = min(min(Q));
U = zeros(n,1); for j=1:n, U(j)=(max(Q(:,j)))-ell;end
% Define the constraints
Constraints = [];
Constraints = [Constraints, sum(x) == 1];  % Sum of x equals 1
Constraints = [Constraints, Q * x <= alpha + z];  % Matrix form for e_j' * Q * x <= alpha + z_j
Constraints = [Constraints, x <= y];                    % Vectorized constraint x_j <= y_j
Constraints = [Constraints, z <= U .* (1 - y)];         % Vectorized constraint z_j <= U_j * (1 - y_j)
Constraints = [Constraints, x >= 0, z >= 0];            % Non-negativity constraints

% Define the objective
Objective = alpha;

% Set up options and solve
options = sdpsettings('verbose', 0, 'solver', 'gurobi', 'gurobi.TimeLimit',3000,'savesolveroutput',1);
sol = optimize(Constraints, Objective, options);
% Check and display results
 x_value = value(x);
 alpha_value = value(alpha);
 Time = sol.solvertime;
 LB= sol.solveroutput.result.objbound;
 UB=value(x)'*Q*value(x);
 disp('____________MILO-StQO__________'); 
 fprintf('UB MILO:   %s\n', mat2str(UB));
 fprintf('LB MILO:   %s\n', mat2str(LB));
 fprintf('Time MILO: %s\n', mat2str(Time));
 disp('_____________________________');   
end
