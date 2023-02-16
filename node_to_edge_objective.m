function [output] = node_to_edge_objective(x,N,prob,xi,alpha)

r = alpha.*x;     % users' reward

% objective function to be maximized
sum = 0;
for i = 1:N
    sum = sum + prob*(x(i) - xi*r(i));
end

% objective function to be minimized
output = -sum;

end

