function [LB,UB,Total_Time,sol_x] = Static_ARO_Concave_QO(Q, c, A, b)
% Code written by Abbas Khademi 
% based on the paper: 
% "Quadratic Optimization Through the Lens of Adjustable Robust Optimization"
% by Abbas Khademi & Ahmadreza Marandi.
% ARO-QO Algorithm: with full static decision rule for Concave Quadratic Minimization
% This function solves a concave quadratic minimization problem with constraints.
% The problem is defined as:
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

% Requirements: Yalmip & gurobi

format long
[m_x,n_x]=size(A);

% Upper Bound (UB) Step: collect scenarios based on robust counterpart
clear('yalmip');
x=sdpvar(n_x,1);
sol=optimize([A*x>=b, x>=0],(0.5*c')*x ,sdpsettings('verbose', 0,'solver','gurobi'));
my_time =  sol.solvertime;
UB= value(x'*Q*x+c'*x);
x0=value(x);                   % Collect scenario
S=[];
Qx=[];
cx=c'*x0;
S=[S,x0];
%
for i=1:n_x
    clear('yalmip');
    x=sdpvar(n_x,1);
    sol=optimize([A*x>=b, x>=0], Q(i,:)*x ,sdpsettings('verbose', 0,'solver','gurobi'));
    my_time = my_time + sol.solvertime;
    x0=value(x);
    Qx=[Qx;Q(i,:)*x0];
    UB2= value(x'*Q*x+c'*x);
    if UB2<UB   % Check the quality of scenarios and select the best one
       UB=UB2;
       S=[S,x0];
    end
end
S_time= my_time ;
fprintf('senario Time:%s\n',mat2str(S_time));
%%
% Mountain Climbing Procedure
x0 =S(:,end);   % Set initial point
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
        fprintf('UB Time: %s\n',mat2str(UB_Time));
        fprintf('UB Value: %s\n',mat2str(UB));
        break
    end
end
    % 
clear('yalmip');
tau=sdpvar;
w=sdpvar(m_x,1);
CC=[0.5*cx + b'*w >= tau, A'*w <= Qx+0.5*c, w>= 0];
p4=optimize(CC,-tau ,sdpsettings('verbose', 0,'solver','gurobi'));
LB=value(tau);
fprintf('LB value: %s\n',mat2str(LB));
Total_Time=p4.solvertime+UB_Time;
fprintf('Total time: %s\n',mat2str(Total_Time));
end
