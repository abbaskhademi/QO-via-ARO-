function [UB, LB, x_can, Time] = Linear_ARO_QO(Q, c, A, b)
% Code written by Abbas Khademi 
% based on the paper: 
% "Quadratic Optimization Through the Lens of Adjustable Robust Optimization"
% by Abbas Khademi & Ahmadreza Marandi.

% Linear ARO-QO; setting u_x=x and w_x=w
% Requirements: Yalmip & Mosek

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
% x_can ........ Candidate solution vector                        %
% Time ......... Total solver time                                %
%=================================================================%
format long
warning off

% Setting parameters
[m, n] = size(A);

%% Collect worst-case scenarios
scenario_value = [];
Time = 0;
Scenario = [];
Qx = [];
for k = 1:n
    clear('yalmip');
    x = sdpvar(n, 1);  
    P1 = optimize([A*x == b, x >= 0], Q(k,:)* x, sdpsettings('verbose',0, 'solver', 'mosek'));
    x = value(x);
    Scenario =[Scenario, x];
    Qx = [Qx; Q(k,:)*x];
    Time = Time + P1.solvertime; % Solver time for finding each scenario 
    scenario_value=[scenario_value,  (x')*Q*x + c'*x];
end
clear('yalmip');
x = sdpvar(n, 1);  
P1 = optimize([A*x == b, x >= 0], 0.5*c'*x, sdpsettings('verbose',0, 'solver', 'mosek'));
x = value(x);
cx = c'*x;
Scenario = [Scenario, x];
Time = Time + P1.solvertime;  % Total time to reach scenarios 
scenario_value=[scenario_value, (x')*Q*x + c'*x];
% 

%% Optimization for lower bound
clear('yalmip');
tau = sdpvar;
w = sdpvar(m, 1);
Cons = [];
Cons = [Cons, 0.5*cx + b'*w >= tau, -A'*w + Qx >= -0.5*c];
sol=optimize(Cons, -tau, sdpsettings('verbose', 0, 'solver', 'mosek'));
Time_LB = sol.solvertime;
Time = sol.solvertime + Time;
LB = value(tau);
  

%% Finding the best scenario as initial point
[UB0, k] = sort(scenario_value);
UB = UB0(1);
x_0 = Scenario(:, k(1)); % Best scenario 
% fprintf('senario ARO: %s\n', mat2str(UB));

%% Finding the Candidate solution & upper bound
[x_can, UB, MCP_time] = Mountain_Climbing_Procedure(Q, A, b, c, x_0);
Time = Time + MCP_time;
%

disp('_______Linear_ARO_QO_________');
fprintf('LB:   %s\n', mat2str(LB));
fprintf('UB:   %s\n', mat2str(UB));
fprintf('Time:  %s\n', mat2str(Time));
disp('_____________________________'); 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x_end, UB_end, MCP_time] = Mountain_Climbing_Procedure(Q, A, b, c, x_0)
% Mountain Climbing Procedure (MCP) based on Bi_QO: Improve the initial solution

%____________Bi_QO Problem Definitio_________________________%
% min  0.5(x' Qplus x+c'x + y' Qplus y+c'y)+ x' Qminus y     %
% s.t.      Ax = b, x >= 0,                                  %
%           Ay = b, y >= 0.                                  %
%____________________________________________________________%

%=========================================================================%
% ____________ INPUT:                                                     %
% Q       - Matrix of quadratic coefficients of the objective function.   %
% A       - Matrix for linear constraints.                                %
% b       - Right-hand side vector for constraints.                       %
% c       - Coefficient vector for the linear part of the objective.      %
% x_0     - Initial point for optimization.                               %
% ____________ OUTPUT:                                                    %
% x_end   - Candidate solution vector at the end of MCP.                  %
% UB_end  - Upper bound value of the objective at x_end.                  %
% MCP_time- Total solver computation time.                                %
%=========================================================================%

%%Representation: Q = Qplus + Qminus, where (Qplus) & (-Qminus) are Positive Definite
[Qplus, Qminus] = Representation_2(Q);
% Requirements: Yalmip & Mosek

% Setting parameters
n = size(Q,2);
% Initializing 
y_matrix = [x_0];
y_value = [x_0'*Q*x_0+c'*x_0];
MCP_time = 0; 
while true
    xe = y_matrix(:, end);
    clear('yalmip');
    y = sdpvar(n, 1); % y is an n-dimensional optimization variable
    % Optimization problem
    sol = optimize([A*y == b, y >= 0], 0.5*(y'*Qplus*y+c'*y) + xe'*Qminus*y, sdpsettings('solver', 'mosek', 'verbose', 0));
    y_matrix = [y_matrix, value(y)];
    y_value = [y_value, value(y)' * Q * value(y)+c'*value(y)];
    MCP_time = MCP_time + sol.solvertime; 
    % Convergence check
    if abs(y_value(end) - y_value(end - 1)) < 0.00001
        x_end = y_matrix(:, end);
        UB_end = y_value(end);
        break;
    end
end
%
end


%_________________________________________________________________________%
function [Qplus, Qminus] = Representation_1(Q)
n = size(Q,2);
I = eye(n);
eps = 0.00001;
Qplus = Q - ((min(eig(Q)) - eps) * I); 
Qminus = (min(eig(Q)) - eps) * I;
end
%_________________________________________________________________________%
function [Qplus, Qminus] = Representation_2(Q)
n = size(Q,2);
I = eye(n);
eps = 0.00001;
Qplus = ((max(eig(Q)) + eps) * I);
Qminus = Q - Qplus;
end
%_________________________________________________________________________%
function [Qplus, Qminus] = Representation_3(Q)
n = size(Q,2);
[V,D] = eig(Q);
V = real(V); D = real(D);
[~,ind] = sort(diag(D), 'descend');
k = sum(diag(D)>=0); 
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
%_________________________________________________________________________%