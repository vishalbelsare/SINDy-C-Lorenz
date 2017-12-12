%% Discrete Duffing Equation with SINDy


clear
%% set up system, store actual solution
alpha = 2.75;
beta = 0.2;


f = @(x) [x(2); -beta * x(1) + alpha * x(2) - x(2)^3];
n = 2;      % dimension
m = 20000;  % number of samples

X = zeros(m,n);     % pre-allocate

X(1,:) = f([.9,.01]);
for i = 1:m-1
    X(i+1,:) = f(X(i,:));
end


% 
% % add noise to the derivative data
% eta = 0.5;      % noise magnitude
% X = X + eta * randn(size(X));


%% SINDy Algorithm:

% library of nonlinear functions 
lib = @(t,a) [1; a(1); a(2); a(1)^2; a(2)^2; a(1)*a(2); ...
        a(1)^3; a(2)^3; a(1)^2*a(2); a(2)^2*a(1); ...
        sin(a(1)); cos(a(1)); sin(a(2)); cos(a(2));]';

% input data into library
for i = 1:m-1
  Theta(i,:) = lib(i,X(i,:));
end

Xi = SINDy(Theta, X(2:m,:));        % coefficient matrix


%% Recover data with Xi 

Y = @(i,a) Xi' * lib(i,a)';

Xs(1,:) = Y(1,[.9,.01]);
for i = 1:m-1
    Xs(i+1,:) = Y(i,Xs(i,:));
end


%% plot actual solution
fontsize = 20;
figure(1)
clf
scatter(X(:,1), X(:,2),2,'k');
xlabel('$x_k$','Interpreter','Latex','Fontsize',fontsize)
ylabel('$y_k$','Interpreter','Latex','Fontsize',fontsize)

%% plot SINDy solution
figure(2)
clf
scatter(Xs(:,1), Xs(:,2),2,'k');
xlabel('$x_k$','Interpreter','Latex','Fontsize',fontsize)
ylabel('$y_k$','Interpreter','Latex','Fontsize',fontsize)



%% compute error
Error = SimulationError(X(:,2), Xs(:,2));