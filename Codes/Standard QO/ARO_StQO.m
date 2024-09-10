function [UB, x_can, Time] = ARO_StQO(Q)
% Code written by Abbas Khademi 
% based on the paper: 
% "Quadratic Optimization Through the Lens of Adjustable Robust Optimization"
% by Abbas Khademi & Ahmadreza Marandi.
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
% x_can ........ Candidate solution vector                        %
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

[x_can, UB, Time] = MCP_StQO_x0(Q,Qplus,Qminus, x_k, Time_sena);


%%

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x_end, u_end, time] = MCP_StQO_x0(Q,Qplus,Qminus, x_0, time_x_0)
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
        x_end = y_matrix(:, end);
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