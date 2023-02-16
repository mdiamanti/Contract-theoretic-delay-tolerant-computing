%% Evaluation results w.r.t. different cost and risk aversion

% variables initialization

a_tot = zeros(length(1:1:10),length(10:1:20));
s_tot = zeros(length(1:1:10),length(10:1:20));
t_tot = zeros(length(1:1:10),length(10:1:20));
w_tot = zeros(length(1:1:10),length(10:1:20));
utility_e_tot = zeros(length(1:1:10),length(10:1:20));
utility_f_tot = zeros(length(1:1:10),length(10:1:20));

a_compl_tot = zeros(length(1:1:10),length(10:1:20));
t_compl_tot = zeros(length(1:1:10),length(10:1:20));
utility_e_compl_tot = zeros(length(1:1:10),length(10:1:20));
utility_f_compl_tot = zeros(length(1:1:10),length(10:1:20));

sigma2 = 0.1;                               % edge server's variance of effort
epsilon = normrnd(0,sigma2,1,1);            % edge server's error between effort and performance

i = 1;
j = 1;
for c = 0.1:0.01:0.2
    for eta = 0.1:0.01:0.2

        eta = eta * 10^4;
        
        % incomplete information
        s = 1/(1+eta*c*sigma2);                     % edge server's variable compensation
        t = -s^2/(2*c) + eta*s^2*sigma2/2;          % edge server's fixed compensation
        a = s/c;                                    % edge server's effort
    
        q = a + epsilon;                            % edge server's performance
        w = t + s*q;                                % edge server's overall compensation
        psi = c*a^2/2;                              % edge server's cost of effort

        utility_e = t + s*a - c*a^2/2 - eta*s^2*sigma2/2;     % edge server's utility 
        utility_f = (1-s)*(s/c)-t;                            % fog server's utility
        
        % complete information
        a_compl = 1/c;                              % edge server's effort
        t_compl = 1/(2*c);                          % edge server's fixed compensation

        utility_e_compl = t_compl - c*a_compl^2/2;  % edge server's utility
        utility_f_compl = a_compl - t_compl;        % fog server's utility
        
        % final results
        a_tot(i,j) = a;
        s_tot(i,j) = s;
        t_tot(i,j) = t;
        w_tot(i,j) = w;
        utility_e_tot(i,j) = utility_e;
        utility_f_tot(i,j) = utility_f;

        a_compl_tot(i,j) = a_compl;
        t_compl_tot(i,j) = t_compl;
        utility_e_compl_tot(i,j) = utility_e_compl;
        utility_f_compl_tot(i,j) = utility_f_compl;
        
        j = j+1;
    end   
    i = i+1;
    j = 1;
end

save('e2f_contract_res1.mat','a_tot','s_tot','t_tot','w_tot','utility_f_tot',...
    'a_compl_tot','t_compl_tot','utility_f_compl_tot');

% plots
figure();
surf(a_tot);
xlim([0 10]);
xticks([0 2 4 6 8 10]);
xlabel('Cost of effort');
xticklabels({'0.1','0.12','0.14','0.16','0.18','0.2'});
ylim([0 10]);
yticks([0 2 4 6 8 10]);
yticklabels({'1000','1200','1400','1600','1800','2000'});
ylabel('Risk aversion');
zlabel('Edge server''s effort');
set(gca,'FontSize',14);
grid on;
[caz,cel] = view
v = [-5 -2 5];
[caz,cel] = view(v)

figure();
surf(t_tot);
xlim([1 10]);
xlabel('Edge server''s cost of effort');
xticklabels({'0.12','0.14','0.16','0.18','0.2'});
ylim([1 10]);
yticklabels({'1200','1400','1600','1800','2000'});
ylabel('Edge server''s risk aversion');
zlabel('Edge server''s fixed reward');
set(gca,'FontSize',16);
grid on;

figure();
surf(s_tot);
xlim([1 10]);
xlabel('Edge server''s cost of effort');
xticklabels({'0.12','0.14','0.16','0.18','0.2'});
ylim([1 10]);
yticklabels({'1200','1400','1600','1800','2000'});
ylabel('Edge server''s risk aversion');
zlabel('Edge server''s variable reward');
set(gca,'FontSize',16);
grid on;

figure();
surf(t_tot);
hold on;
surf(s_tot);
xlim([0 10]);
xticks([0 2 4 6 8 10]);
xlabel('Cost of effort');
xticklabels({'0.1','0.12','0.14','0.16','0.18','0.2'});
ylim([0 10]);
yticks([0 2 4 6 8 10]);
yticklabels({'1000','1200','1400','1600','1800','2000'});
ylabel('Risk aversion');
zlabel('Edge server''s reward');
lgd = legend('Fixed reward','Variable reward');
set(lgd,'FontSize',14);
set(gca,'FontSize',14);
set(gca,'zscale','log')
grid on;
[caz,cel] = view
v = [-5 -2 5];
[caz,cel] = view(v)

figure();
surf(utility_f_tot);
xlim([0 10]);
xticks([0 2 4 6 8 10]);
xlabel('Cost of effort');
xticklabels({'0.1','0.12','0.14','0.16','0.18','0.2'});
ylim([0 10]);
yticks([0 2 4 6 8 10]);
yticklabels({'1000','1200','1400','1600','1800','2000'});
ylabel('Risk aversion');
zlabel('Fog server''s utility');
set(gca,'FontSize',14);
grid on;
[caz,cel] = view
v = [-5 -2 5];
[caz,cel] = view(v)
