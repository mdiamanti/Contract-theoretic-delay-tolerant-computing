%% Initializations

N = 10;                         % number of users
%T = randi([50 2000],1,N);      % users' service time requirement [msec]
T = randi([500 2000],1,N);      % users' service time requirement [msec]
B = randi([1000 5000],1,N);     % users' total amount of data [KB]
phi = randi([1000 2000],1,N);   % users' task's intensity [CPU cycles/bytes]

C = phi.*B;                     % users' task's intensity [CPU cycles]

%% User-to-edge server contract - Adverse selection problem

prob = 1/N;                     % users' probability

w1 = 0.6;
w2 = 1-w1;
theta = w1*T/sum(T) + w2*phi/sum(phi);    % users' type
theta = sort(theta);  

kappa = 0.3;                    % users' cost of effort
xi = 0.25;                      % edge server's cost of effort
alpha = 8*sqrt(theta);          % percentage from the edge server's monetary savings

fun = @(x)node_to_edge_objective(x,N,prob,xi,alpha);            % objective function
nonlcon = @(x)node_to_edge_constraint(x,N,kappa,theta,alpha);   % non-linear constraints

lb = zeros(1,N);                % solution's lower bounds
ub = ones(1,N);                 % solution's upper bounds

x0 = 0.1:0.1:1;

options = optimoptions('fmincon','Display','iter','Algorithm','sqp-legacy');    % optimization options

[x,fval,exitflag] = fmincon(fun,x0,[],[],[],[],lb,ub,nonlcon,options);          % call fmincon

p = x;                                      % users' performance
r = alpha.*x;                               % users' reward 

utility_n = theta.*sqrt(r) - kappa*p;       % users' utility
utility_e = prob * (p-xi*r);                % edge servers' utility

% final metrics
c = p.*C;                                   % users' amount of task to be offloaded at the fog
b = c./phi;                                 % users' amount of data to be offloaded at the fog

D = sum(b);                                 % edge server's total amount of data that can be potentially offloaded at the fog

%% Adverse selection - Complete information contract

p_compl = theta.^2/(2*xi*kappa^2);          % users' performance
r_compl = (theta/(2*xi*kappa)).^2;          % users' reward

% final metrics
c_compl = p_compl.*C;                       % users' amount of task to be offloaded at the fog
b_compl = c_compl./phi;                     % users' amount of data to be offloaded at the fog

D_compl = sum(b_compl);                     % edge server's total amount of data that can be potentially offloaded at the fog

%% Adverse selection - Uniform linear-pricing contract

prob = 1/N;                     % users' probability

w1 = 0.6;
w2 = 1-w1;
theta = w1*T/sum(T) + w2*phi/sum(phi);    % users' type
theta = sort(theta);  

kappa = 0.3;                    % users' cost of effort
xi = 0.25;                      % edge server's cost of effort
alpha = 8*sqrt(sum(theta)/N);   % percentage from the edge server's monetary savings

fun = @(x)node_to_edge_objective(x,N,prob,xi,alpha);            % objective function
nonlcon = @(x)node_to_edge_constraint(x,N,kappa,theta,alpha);   % non-linear constraints

lb = zeros(1,N);                % solution's lower bounds
ub = ones(1,N);                 % solution's upper bounds

x0 = 0.1:0.1:1;

options = optimoptions('fmincon','Display','iter','Algorithm','sqp-legacy');    % optimization options

[x,fval,exitflag] = fmincon(fun,x0,[],[],[],[],lb,ub,nonlcon,options);          % call fmincon

p_lin = x;                                      % users' performance
r_lin = alpha.*x;                               % users' reward 

utility_n_lin = theta.*sqrt(r_lin) - kappa*p_lin;   % users' utility
utility_e_lin = prob * (p_lin-xi*r_lin);            % edge servers' utility

% final metrics
c_lin = p_lin.*C;                               % users' amount of task to be offloaded at the fog
b_lin = c_lin./phi;                             % users' amount of data to be offloaded at the fog

D_lin = sum(b_lin);                             % edge server's total amount of data that can be potentially offloaded at the fog

%% Edge-to-fog server contract - Moral hazard problem

sigma2 = 0.1;                               % edge server's variance of effort
epsilon = normrnd(0,sigma2,1,1);            % edge server's error between effort and performance

c = 10;                                     % edge server's unit cost of effort [0.1,0.2]
eta = 10;                                   % edge server's coefficient of CARA [1000,2000]

s = 1/(1+eta*c*sigma2);                     % edge server's variable compensation
t = -s^2/(2*c) + eta*s^2*sigma2/2;          % edge server's fixed compensation
a = s/c;                                    % edge server's effort

q = a + epsilon;                            % edge server's performance
w = t + s*q;                                % edge server's overall compensation
psi = c*a^2/2;                              % edge server's cost of effort

utility_e = t + s*a - c*a^2/2 - eta*s^2*sigma2/2;     % edge server's utility 
utility_f = (1-s)*(s/c)-t;                            % fog server's utility

%-----------------------other formulation----------------------------------
c = 10;                                     % edge server's unit cost of effort [0.1,0.2]
eta = 15;                                   % edge server's coefficient of CARA [1000,2000]

mu = 20;
s = 1/(1+mu*eta*c*sigma2);                     % edge server's variable compensation
t = -mu^2*s^2/(2*c) + eta*mu^2*s^2*sigma2/2;   % edge server's fixed compensation
a = mu*s/c;                                    % edge server's effort

q = a + epsilon;                            % edge server's performance
w = t + s*q;                                % edge server's overall compensation
psi = c*a^2/2;                              % edge server's cost of effort

utility_e = t + mu*s*a - c*a^2/2 - eta*mu^2*s^2*sigma2/2;     % edge server's utility 
utility_f = (1-mu*s)*(mu*s/c)-t;                            % fog server's utility

%% Moral hazard - Complete information contract

a_compl = 1/c;                              % edge server's effort
t_compl = 1/(2*c);                          % edge server's fixed compensation

utility_e_compl = t - c*a^2/2;              % edge server's utility
utility_f_compl = a_compl - t_compl;        % fog server's utility


