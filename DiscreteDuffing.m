clear


alpha = 2.75;
beta = 0.2;

f = @(x) [x(2); -beta * x(1) + alpha * x(2) - x(2)^3];

n = 2;
m = 20000;

X = zeros(m,n);

X(1,:) = f([.9,.01]);
for i = 1:m-1
    X(i+1,:) = f(X(i,:));
end

% disp(X(:,1));

figure()
clf
scatter(X(:,1), X(:,2),1);






%% SINDy Algorithm:

% library of nonlinear functions 
lib = @(t,a) [1; a(1); a(2); a(1)^2; a(2)^2; a(1)*a(2); ...
        a(1)^3; a(2)^3; a(1)^2*a(2); a(2)^2*a(1); ...
        sin(a(1)); cos(a(1)); sin(a(2)); cos(a(2));]';

% input data into library
for i = 1:m-1
  Theta(i,:) = lib(i,X(i,:));
end

Xi = SINDy(Theta, X(2:m,:));

% now have solution Xi that describes which nonlinear terms are active in
% the dynamics of the input system


%% Recover data with Xi 

Y = @(i,a) Xi' * lib(i,a)';

Xs(1,:) = Y(1,[.9,.01]);
for i = 1:m-1
    Xs(i+1,:) = Y(i,Xs(i,:));
end

figure()
clf
scatter(Xs(:,1), Xs(:,2),2);

Error = SimulationError(X(:,2), Xs(:,2));