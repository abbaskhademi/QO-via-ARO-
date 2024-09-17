% Example data generation for problem #71
% Selvi et al. (2022) consider the following problem and data generation:
%================================%
%      max    x'Px + ell'x       %
%      s.t.   Dx <= d,           %
%             x >= 0.            %
%================================%
m = 800;  % Number of constraints
n = 700;  % Dimension of the variable x
rng(23);  % Set seed for the random number generator to a specific value
ell = zeros(n,1);               % Objective is x^T L^T L x + ell^T x
L = rand(n,n);                  % Generate L
P = L'*L;                       % P is L^T L 
D = eye(n);                     % First half of constraints are x <= x_u
D = [D; rand(m-n,n)*2];         % Second half of the constraints are Dx <= d
d = 8*ones(n,1);                % x_u = 8 
d = [d; randi([30,60],m-n,1)];  % RHS of Dx <= d
% Changing the sign of data leads to a concave minimization problem as below.
Q = -P; A = -D; b = -d; c = -ell;
%================================%
%      min    x'Qx + c'x         %
%      s.t.   Ax >= b,           %
%             x >= 0.            %
%================================%
save('#71 Problem_mx800nx700.mat','Q','A','b','c');
