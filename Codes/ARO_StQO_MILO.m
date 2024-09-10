%function [LB,UB,Time]=ARO_cut_MILO2_cplex(Q,UB,T0)
function [LB,UB,Time]=ARO_StQO_MILO(Q)

[UB, T0] = ARO_StQO(Q);
n = size(Q, 1); 
M = zeros(n,1); for j=1:n, M(j)=(max(Q(:,j)))-min(min(Q));end
%%    find warm start
yalmip clear; 
% Define the decision variables
x = sdpvar(n, 1);          % Vector x
z = sdpvar(n, 1);          % Vector z
alpha = sdpvar(1);         % Scalar alpha
y = binvar(n, 1); 
Constraints = [Q * x <= alpha + z, sum(x) == 1, 0 <= x <= y, 0 <= z <= M .* (1 - y),  alpha <= UB ];
Sol = optimize(Constraints, 0, sdpsettings('solver', 'cplex', 'verbose', 0)); 
T1=Sol.solvertime;
fprintf('Time finding warm: %s\n', mat2str(T1));
disp(yalmiperror(Sol.problem))
aa =value(alpha);
sum(value(x));
y_x=value(y);
z_x=value(z);
xx=value(x);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         solve
yalmip clear; 
x = sdpvar(n, 1);          % Vector x
z = sdpvar(n, 1);          % Vector z
y = binvar(n, 1);          % Binary vector y
alpha = sdpvar(1);         % Scalar alpha
% Define the constraints
Constraints = [];
Constraints = [Constraints, sum(x) == 1];  % Sum of x equals 1
Constraints = [Constraints, Q * x <= alpha + z];  % Matrix form for e_j' * Q * x <= alpha + z_j
Constraints = [Constraints, x <= y];                    % Vectorized constraint x_j <= y_j
Constraints = [Constraints, z <= M .* (1 - y)];         % Vectorized constraint z_j <= U_j * (1 - y_j)
Constraints = [Constraints, x >= 0, z >= 0];            % Non-negativity constraints
Constraints = [Constraints, min(min(Q))<= alpha];            % Non-negativity constraints

% Define the objective
Objective = alpha;
warmstart(x,xx);warmstart(alpha,aa);warmstart(y,y_x);warmstart(z,z_x);
% Set up options and solve
options = sdpsettings('verbose',0, 'solver', 'gurobi','warmstart',1,'gurobi.TimeLimit',max([3000-T1-T0,1]),'savesolveroutput',1);  % Using Gurobi solver; change as needed
sol = optimize(Constraints, Objective, options);
T2=sol.solvertime
% Display results
disp(yalmiperror(sol.problem))
UB=value(x)'*Q*value(x);
fprintf('UB MILO: %s\n', mat2str(UB));
LB= sol.solveroutput.result.objbound;
fprintf('LB: %s\n', mat2str(LB));
Time = T0+T1+T2;
fprintf('Time MILO: %s\n', mat2str(Time));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%5
function [UB, Time] = ARO_StQO(Q)
format long
warning off

%__________StQO Problem Definition__________%
% min    x'Qx                               %
% s.t.   e'x = 1, x >= 0.                   %
%___________________________________________%

%=================================================================%
% ____________ INPUT:                                             %
% Q ........... Matrix of the objective function                  %
% ____________ OUTPUT:                                            %
% UB ........... Upper bound value                                %
% Time ......... Total solver time                                %
%=================================================================%

% Requirements: Yalmip & Mosek

% Setting parameters
[~, n] = size(Q);
I = eye(n);

%% Representation: Q = Qplus + Qminus, where Qplus & -Qminus are Positive Definite

% [Qplus, Qminus] = Representation_1(Q); R=1; % Representation 1
 [Qplus, Qminus] = Representation_2(Q); R=2; % Representation 2
% [Qplus, Qminus] = Representation_3(Q); R=3; % Representation 3


%% Collect worst-case scenarios
S = zeros(n, 1);
Time = zeros(n, 1);
Scenario = zeros(n, n);
for k = 1:n
    clear('yalmip');
    x = sdpvar(n, 1);  
    q = Qminus(k, :);
    P1 = optimize([sum(x) == 1, x >= 0], 0.5 * (x') * Qplus * x + q * x, sdpsettings('verbose', 0, 'solver', 'mosek'));
    %
    x = value(x);
    Scenario(:, k) = x;
    Time(k) = P1.solvertime; % Solver time for each scenario
    S(k) = 0.5 * (x') * Qplus * x + q * x;
end
Time_sena = sum(Time); 
% Add vertieces
Scenario = [Scenario,I]; 

% Finding the best scenario and upper bound
upper = zeros(max(size(Scenario)), 1);
for i = 1:max(size(Scenario)) 
    upper(i) = (Scenario(:, i))' * Q * (Scenario(:, i));
end
[UB0, k] = sort(upper);
UB = UB0(1);
x_k = Scenario(:, k(1)); % Best scenario L2_StQO

% Output display
fprintf('UB Best Scenario: %s\n', mat2str(UB));
fprintf('Scenario Time: %s\n', mat2str(Time_sena));
%%

[UB, Time] = MCP_StQO_x0(Q,Qplus,Qminus, x_k, Time_sena);


%%

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [u_end, time] = MCP_StQO_x0(Q,Qplus,Qminus, x_0, time_x_0)
% Mountain Climbing Procedure (MCP) based on Bi_StQO: Improve the initial solution
format long
warning off

%________________Bi_StQP Problem Definition_________________%
% min   0.5 * {x' * (Q+) * x + y' * (Q+) * y} + x' * (Q-) * y %
% s.t.   e' * x = 1, x >= 0, e' * y = 1, y >= 0.               %
%___________________________________________________________%

%====================================================%
% ____________ INPUT:                                %
% Q .......... Matrix of the objective function      %
% x_0 ......... Initial point                        %
% time_x_0 .... Time to reach initial point          %
% ____________ OUTPUT:                               %
% x_end ....... Candidate solution                   %
% u_end ....... Upper bound value                    %
% time ........ Total solver times to reach UB       %
%====================================================%

% Requirements: Yalmip & Mosek

% Setting parameters
[~, n] = size(Q);
I = eye(n);


% Initializing variables for MCP
y_matrix = [x_0];
y_value = [x_0' * Q * x_0];
time_y_cal = 0; 

% MCP iteration
while true
    xe = y_matrix(:, end);
    clear('yalmip');
    y = sdpvar(n, 1); % y is an n-dimensional optimization variable

    % Optimization problem
    sol = optimize([sum(y) == 1, y >= 0], 0.5 * y' * Qplus * y + xe' * Qminus * y, sdpsettings('solver', 'mosek', 'verbose', 0));
    y_matrix = [y_matrix, value(y)];
    y_value = [y_value, value(y)' * Q * value(y)];
    time_y_cal = time_y_cal + sol.solvertime;

    % Convergence check
    if abs(y_value(end) - y_value(end - 1)) < 0.000001
        %x_end = y_matrix(:, end);
        u_end = y_value(end);
        time = time_y_cal + time_x_0; % Total time
        %%
        fprintf('UB MC: %s\n', mat2str(u_end));
        fprintf('Total Time: %s\n', mat2str(time));
        %%
        break;
    end
end
end
%%%%%%%%%%%%%%%%%%%%55
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Qplus, Qminus] = Representation_3(Q)
n = size(Q,2);
[V,D] = eig(Q);
V = real(V); D = real(D);
[~,ind] = sort(diag(D), 'descend');
k = sum(diag(D)>=0); %number of nonzero elements
D = D(ind,ind);
V  = V(:,ind);
eps = 0.00001;
D_plus_diag = diag(D);
D_plus_diag(k+1:end) = 0;
D_plus = diag(D_plus_diag);
D_plus = D_plus + (eye(n)*eps);
%
D_minus_diag = diag(D);
D_minus_diag(1:k) = 0;
D_minus = diag(-D_minus_diag);
D_minus = D_minus+(eye(n)*eps);
Qplus = V*D_plus*V';
Qminus = -V*D_minus*V';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Qplus, Qminus] = Representation_2(Q)
n = size(Q,2);
I = eye(n);
Qplus = ((max(eig(Q)) + 0.00001) * I);
Qminus = Q - Qplus;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Qplus, Qminus] = Representation_1(Q)
n = size(Q,2);
I = eye(n);
Qplus = Q - ((min(eig(Q)) - 0.00001) * I); 
Qminus = (min(eig(Q)) - 0.00001) * I;
end