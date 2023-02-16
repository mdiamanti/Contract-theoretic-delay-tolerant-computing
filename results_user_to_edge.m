%% Evaluation results w.r.t. user types

% variables initialization
theta_avg = 0;

p_avg = 0;
r_avg = 0;

utility_n_avg = 0;
utility_e_avg = 0;
utility_n_test_avg = 0;

D_avg = 0;

p_compl_avg = 0;
r_compl_avg = 0;

utility_n_compl_avg = 0;
utility_e_compl_avg = 0;

D_compl_avg = 0;

p_lin_avg = 0;
r_lin_avg = 0;

utility_n_lin_avg = 0;
utility_e_lin_avg = 0;

for ite = 1:50
    
    N = 10;                         % number of users
    T = randi([500 2000],1,N);      % users' service time requirement [msec]
    B = randi([1000 5000],1,N);     % users' total amount of data [KB]
    phi = randi([1000 2000],1,N);   % users' task's intensity [CPU cycles/bytes]

    C = phi.*B;                     % users' task's intensity [CPU cycles]

    prob = 1/N;                     % users' probability
    
    w1 = 0.6;
    w2 = 1-w1;
    theta = w1*T/sum(T) + w2*phi/sum(phi);    % users' type
    theta = sort(theta);               

    kappa = 0.3;                    % users' cost of effort
    xi = 0.25;                      % edge server's cost of effort
    alpha = 8*sqrt(theta);          % percentage from the edge server's monetary savings
    
    % incomplete information
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
    
    for i = 1:N
        for j = 1:N
            utility_n_test(i,j) = theta(i)*sqrt(r(j)) - kappa*p(j);
        end 
    end
    
    % final metrics
    c = p.*C;                                   % users' amount of task to be offloaded at the fog
    b = c./phi;                                 % users' amount of data to be offloaded at the fog

    D = sum(b);                                 % edge server's total amount of data that can be potentially offloaded at the fog

    % complete information
    p_compl = theta.^2/(2*xi*kappa^2);          % users' performance
    r_compl = (theta/(2*xi*kappa)).^2;          % users' reward
    
    utility_n_compl = theta.*sqrt(r_compl) - kappa*p_compl;      % users' utility
    utility_e_compl = prob * (p_compl-xi*r_compl);               % edge servers' utility
    
    % final metrics
    c_compl = p_compl.*C;                        % users' amount of task to be offloaded at the fog
    b_compl = c_compl./phi;                      % users' amount of data to be offloaded at the fog

    D_compl = sum(b_compl);                      % edge server's total amount of data that can be potentially offloaded at the fog

    % uniform linear pricing
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
    
    % sum results for averaging
    theta_avg = theta_avg + theta;
    
    p_avg = p_avg + p;
    r_avg = r_avg + r;
    
    utility_n_avg = utility_n_avg + utility_n;
    utility_e_avg = utility_e_avg + utility_e;
    utility_n_test_avg = utility_n_test_avg + utility_n_test;
    
    p_compl_avg = p_compl_avg + p_compl;
    r_compl_avg = r_compl_avg + r_compl;
    
    utility_n_compl_avg = utility_n_compl_avg + utility_n_compl;
    utility_e_compl_avg = utility_e_compl_avg + utility_e_compl;
    
    p_lin_avg = p_lin_avg + p_lin;
    r_lin_avg = r_lin_avg + r_lin;
    
    utility_n_lin_avg = utility_n_lin_avg + utility_n_lin;
    utility_e_lin_avg = utility_e_lin_avg + utility_e_lin;
    
    D_avg = D_avg + D;
    D_compl_avg = D_compl_avg + D_compl;
    
end

theta_avg = theta_avg/ite;

p_avg = p_avg/ite;
r_avg = r_avg/ite;

utility_n_avg = utility_n_avg/ite;
utility_e_avg = utility_e_avg/ite;
utility_n_test_avg = utility_n_test_avg/ite;

p_compl_avg = p_compl_avg/ite;
r_compl_avg = r_compl_avg/ite;

utility_n_compl_avg = utility_n_compl_avg/ite;
utility_e_compl_avg = utility_e_compl_avg/ite;

p_lin_avg = p_lin_avg/ite;
r_lin_avg = r_lin_avg/ite;

utility_n_lin_avg = utility_n_lin_avg/ite;
utility_e_lin_avg = utility_e_lin_avg/ite;

D_avg = D_avg/ite;
D_compl_avg = D_compl_avg/ite;

save('n2e_contract_res1.mat','N','theta_avg','p_avg','r_avg','utility_n_avg',...
    'utility_e_avg','utility_n_test_avg','p_compl_avg','r_compl_avg',...
    'utility_n_compl_avg','utility_e_compl_avg','p_lin_avg','r_lin_avg',...
    'utility_n_lin_avg','utility_e_lin_avg');

% plots
figure();
plot(p_avg,'-d', 'LineWidth',1);
hold on;
plot(p_compl_avg,'-*','LineWidth',1);
hold on;
plot(p_lin_avg,'-o','LineWidth',1);
xlim([1 N]);
xlabel('User''s index','FontSize',18);
ylabel('User''s effort','FontSize',18);
lgd = legend('Incomplete information','Complete information','Uniform contract');
set(lgd,'FontSize',14);
set(gca,'FontSize',18);
grid on;

figure();
plot(r_avg,'-d', 'LineWidth',1);
hold on;
plot(r_compl_avg,'-*','LineWidth',1);
hold on;
plot(r_lin_avg,'-o','LineWidth',1);
xlim([1 N]);
xlabel('User''s index','FontSize',18);
ylabel('User''s reward','FontSize',18);
lgd = legend('Incomplete information','Complete information','Uniform contract');
set(lgd,'FontSize',14);
set(gca,'FontSize',18);
grid on;

utility_n_compl_avg = zeros(1,N);
figure();
plot(utility_n_avg,'-d', 'LineWidth',1);
hold on;
plot(utility_n_compl_avg,'-*','LineWidth',1);
hold on;
plot(utility_n_lin_avg,'-o','LineWidth',1);
xlim([1 N]);
xlabel('User''s index','FontSize',18);
ylabel('User''s utility','FontSize',18);
lgd = legend('Incomplete information','Complete information','Uniform contract');
set(lgd,'FontSize',14);
set(gca,'FontSize',18);
grid on;

figure();
plot(utility_e_avg,'-d', 'LineWidth',1);
hold on;
plot(utility_e_compl_avg,'-*','LineWidth',1);
hold on;
plot(utility_e_lin_avg,'-o','LineWidth',1);
xlim([1 N]);
xlabel('User''s index','FontSize',18);
ylabel('Edge server''s utility','FontSize',18);
lgd = legend('Incomplete information','Complete information','Uniform contract');
set(lgd,'FontSize',14);
set(gca,'FontSize',18);
grid on;

figure();
plot(utility_n_avg+utility_e_avg,'-d', 'LineWidth',1);
hold on;
plot(utility_n_compl_avg+utility_e_compl_avg,'-*','LineWidth',1);
hold on;
plot(utility_n_lin_avg+utility_e_lin_avg,'-o','LineWidth',1);
xlim([1 N]);
xlabel('User''s index','FontSize',20);
ylabel('Social welfare','FontSize',20);
lgd = legend('Incomplete information','Complete information','Uniform contract');
set(lgd,'FontSize',14);
set(gca,'FontSize',20);
grid on;

figure(); 
plot(1:N,utility_n_test_avg(3,1:N),'-d', 'LineWidth',1);
hold on;
plot(1:N,utility_n_test_avg(5,1:N),'-*', 'LineWidth',1);
hold on;
plot(1:N,utility_n_test_avg(7,1:N),'-s', 'LineWidth',1);
xlim([1 N]);
xlabel('User''s index','FontSize',18);
ylabel('User''s utility','FontSize',18);
set(gca,'FontSize',18);
YL = get(gca, 'ylim');
YL(2) = utility_n_test_avg(5,5);
line([5, 5], YL, 'Color',[0.8500 0.3250 0.0980],'LineStyle','--','LineWidth',1);
YL = get(gca, 'ylim');
YL(2) = utility_n_test_avg(3,3);
line([3, 3], YL, 'Color',[0 0.4470 0.7410],'LineStyle','--','LineWidth',1);
YL = get(gca, 'ylim');
YL(2) = utility_n_test_avg(7,7);
line([7, 7], YL, 'Color',[0.9290 0.6940 0.1250],'LineStyle','--','LineWidth',1);
lgd = legend('User 3','User 5','User 7');
set(lgd,'FontSize',14);
grid on;
