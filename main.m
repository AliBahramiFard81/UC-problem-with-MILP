%% by Ali Bahrami Fard
% alibt1313@gmail.com
% UNIT COMMMITMENT
clc;
clear;

%% Inputs
syms p1 p2 p3 L; 

H = [510 + 7.2 * p1 + 0.00142 * p1 ^ 2;
     310 + 7.85 * p2 + 0.00194 * p2 ^ 2;
     78 + 7.97 * p3 + 0.00482 * p3 ^ 2];

H_cost = [0.9;
           1;
           1;];

power_limits = [150,600;
                100,400;
                50,200;];

start_up_cost = [500;500;500];

init_state = [2;-8;-8];

min_up_down_time = [2,2;
                    2,2;
                    2,2;];

demands = [650,1050,700,1100,350,300,980,800,660,310];

vars = symvar(H);

%% Calculating Functions
F = H;
for i=1:length(H)
   F(i,1) = H(i,1) * H_cost(i,1);
end

%% Piecewise linear const function
number_of_pieces = 11;
pieced_power_limit = [];
index = 1;
for i=1:length(power_limits)
    power_limit = power_limits(i,:);
    power_difference = (power_limit(1,2) - power_limit(1,1)) / number_of_pieces;
    prev_differencce = power_limit(1,1);
    for x=1:number_of_pieces + 1
        pieced_power_limit(index,x) = prev_differencce;
        prev_differencce = prev_differencce + power_difference;
    end
    index = index + 1;
end

%% Calculating F
F_piece = [];
for i=1:length(pieced_power_limit(:,1))
    for x=1:length(pieced_power_limit(1,:))
       F_piece(i,x) = subs(F(i,1),pieced_power_limit(i,x));
    end
end

%% Calculating S
S = [];
for i=1:length(pieced_power_limit(:,1))
    power_limit = power_limits(i,:);
    power_difference = (power_limit(1,2) - power_limit(1,1)) / number_of_pieces;
    for x=1:length(pieced_power_limit(1,:)) - 1
        S(i,x) = (F_piece(i,x + 1) - F_piece(i,x)) / power_difference;
    end
end
disp(pieced_power_limit)
disp(S)

%% Calculating F for each piece
F_x = [];
for i=1:length(pieced_power_limit(:,1))
    for x=1:length(pieced_power_limit(i,:))
        F_x(i,x) = double(subs(F(i,1),pieced_power_limit(i,x)));
    end
end
F_min = F_x(:,end);
disp(F_x)
disp(F_min)

%% Solving MILP

% number of thermal units
N = length(H);
% number of time steps
T = length(demands);

U = optimvar('U',N,T,'Type','integer','LowerBound',0,'UpperBound',1);
Y = optimvar('Y',N,T,'Type','integer','LowerBound',0,'UpperBound',1);
Z = optimvar('Z',N,T,'Type','integer','LowerBound',0,'UpperBound',1);
Pit = optimvar('Pit',N,T);

problem = optimproblem("Objective", sum(F_min'*U) + sum(start_up_cost'*Y));

% unit generation to load constrain
problem.Constraints.pd = sum(Pit,1) == demands;

% unit power min max constrains
problem.Constraints.min = Pit >= repmat(power_limits(:,1),1,T).*U; 
problem.Constraints.max = Pit <= repmat(power_limits(:,2),1,T).*U; 

constrains_UYZ = optimconstr(N, T);
for i=2:T
    constrains_UYZ(:,i-1) = U(:,i) - U(:,i-1) == Y(:,i) - Z(:,i);
end
problem.Constraints.UYZ = constrains_UYZ;

problem.Constraints.YZ = Y + Z <= 1;

% min up/down time constrains

% Minimum Up-Time Constraints
for i = 1:N
    for t = 1:(T - min_up_down_time(i,1) + 1) % Ensure we don't exceed time horizon
        problem.Constraints.(['min_up_time_' num2str(i) '_' num2str(t)]) = ...
            sum(U(i, t:t+min_up_down_time(i,1)-1)) >= min_up_down_time(i,1) * Y(i,t);
    end
end

% Minimum Down-Time Constraints
for i = 1:N
    for t = 1:(T - min_up_down_time(i,2) + 1) % Ensure we don't exceed time horizon
        problem.Constraints.(['min_down_time_' num2str(i) '_' num2str(t)]) = ...
            sum(1 - U(i, t:t+min_up_down_time(i,2)-1)) >= min_up_down_time(i,2) * Z(i,t);
    end
end


options = optimoptions('intlinprog');
[sol,TotalCost,exitflag,output] = solve(problem, 'Options', options);
F_op = [];
for i=1:length(round(sol.Pit(:,1)))
    for x=1:length(round(sol.Pit(i,:)))
        F_op(i,x) = double(subs(F(i,1),round(sol.Pit(i,x))));
    end
end
disp(round(sol.U));
disp(round(sol.Y));
disp(round(sol.Z));
disp(round(sol.Pit));
disp(F_op)