function [c,ceq] = node_to_edge_constraint(x,N,kappa,theta,alpha)

r = alpha.*x;     % users' reward

% non-linear equality constraints
ceq = zeros(1,N);

ceq(1) = theta(1)*sqrt(r(1)) - kappa*x(1);

for i = 2:N
    ceq(i) = theta(i)*sqrt(r(i)) - kappa*x(i) - theta(i)*sqrt(r(i-1))...
        + kappa*x(i-1);
end

% non-linear inequality constraints
c = zeros(1,N);

c(1) = -r(1);

for i = 2:N
    c(i) = r(i-1)-r(i);
end

end

