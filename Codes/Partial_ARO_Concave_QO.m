function [LB,UB,Total_Time,sol_x]=Partial_ARO_Concave_QO(Q,c,A,b) 

% ARO-QO Algorithm: with partial affine decision rule for Concave Quadratic Mimimization

% This function solves the following problem:
%===============================%
%	    min    x'Qx + c'x       %
%	    s.t.     Ax >= b,       %
%                 x >= 0.       %
%===============================%

%=================================================================%
% ____________ INPUT:                                             %
% Q ......... the matrix associated with the objective function   %
% c ......... a column vector associated with the linear part of  %
% the objective function                                          %
% A ......... the matrix associated with the constraints          %
% b .........  a column vector associated RHS of the constraints  %
%                                                                 %
% ____________ OUTPUT:                                            %
% LB ........ Lower Bound                                         %
% UB ........ Upper Bound                                         %
% sol_x ..... Candidate solution based ARO QO                     %
% Total_Time ...... Total time to reach to Lower and Upper bounds %
%=================================================================%


% Requirements: Yalmip & Gurobi

format long
[m_x,n_x]=size(A);
%%
clear('yalmip');
% Define variables of LB problem
tau = sdpvar;
z = sdpvar(m_x,1);           
Z = sdpvar(m_x,n_x,'full');  
alpha = sdpvar(m_x,1); 
Pi = sdpvar(m_x,n_x,'full');  
T =  sdpvar(m_x,m_x,'full');
%%
el = round(m_x/7)+1;
% Define constraints
C = [];
C = C + [Z(el:m,:)==0]; % hybrid
C = C + [b'*(alpha+z)>=tau];
C = C + [A'*alpha<=0.5*c+Z'*b];
C = C + [(b'*T)+z'>=0];
C = C + [A'*T<=Z'];
C = C + [A'*Pi<=(Q-A'*Z)'];
C = C + [b'*Pi+(0.5*c')-(A'*z)'>=0];
C = C + [alpha>=0, Pi(:)>=0, T(:)>=0];
%%
% Solve the LB problem
sol = optimize(C,-tau,sdpsettings('verbose', 0,'solver','gurobi')); 
% Analyze error flags
if sol.problem ~= 0
    disp("Error in optimization:");
    yalmiperror(sol.problem)
else 
% Extract LB value 
    LB_Time = sol.solvertime;  
    LB = value(tau);     
    Algorithm = {'(ARO QO): LB with Partial Decision Rule'};
    On_QP = table(Algorithm,LB,LB_Time)
    Z = value(Z); tau = value(tau); z = value(z);    % the UB problem optimal variables
end


%%
% UB Step: collect scenarios based on robust conterpart
clear('yalmip');
x=sdpvar(n_x,1);
sol=optimize([A*x>=b, x>=0],(0.5*c'+b'*Z)*x ,sdpsettings('verbose', 0,'solver','gurobi'));
my_time =  sol.solvertime;
UB= value(x'*Q*x+c'*x);
x0=value(x);                   % collect scenario
S=[];
S=[S,x0];
%
for i=1:n_x
    clear('yalmip');
    x=sdpvar(n_x,1);
    B= Q-A'*Z;
    sol=optimize([A*x>=b, x>=0], B(i,:)*x ,sdpsettings('verbose', 0,'solver','gurobi'));
    my_time = my_time + sol.solvertime;
    x0=value(x);
    UB2= value(x'*Q*x+c'*x);
    if UB2<UB   % check the quality of senarios and selected best one
       UB=UB2;
       S=[S,x0];
    end
end
%
for i=1:el
    clear('yalmip');
    x=sdpvar(n_x,1);
    sol=optimize([A*x>=b, x>=0], Z(i,:)*x ,sdpsettings('verbose', 0,'solver','gurobi'));
    my_time = my_time + sol.solvertime;
    x0=value(x);
    UB2= value(x'*Q*x+c'*x);
    if UB2<UB   % check the quality of senarios and selected best one
       UB=UB2;
       S=[S,x0];
    end
end
%
%%
% Mountain Climbing Procedure
x0 =S(:,end);   %  Set initial point
while 1
    clear('yalmip');
    y=sdpvar(n_x,1);
    p3=optimize([A*y>=b, y>=0],x0'*Q*y+0.5*c'*y ,sdpsettings('verbose', 0,'solver','gurobi'));
    my_time = my_time + p3.solvertime;
    x0=value(y);
    S=[S,x0];
    if abs(S(:,end)'*Q*S(:,end)+c'*S(:,end)-S(:,end-1)'*Q*S(:,end-1)-c'*S(:,end-1))<=0.00001
        UB= (S(:,end)'*Q*S(:,end)+c'*S(:,end));
        sol_x = S(:,end);    % Solution candidate
        UB_Time=my_time;
        Algorithm = {'(ARO_QO): UB'};
        On_QP = table(Algorithm,UB,UB_Time)
        Total_Time = my_time + LB_Time 
        S;
        break
    end
end
    %            
end
