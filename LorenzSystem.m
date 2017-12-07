clear

%% Set up and solve Lorenz system to get state X and dX

sigma=10;
beta=8/3;
rho=28;
tmax = 80;

% Lorenz equations
f = @(t,a) [-sigma*a(1) + sigma*a(2); rho*a(1) - a(2) - a(1)*a(3); -beta*a(3) + a(1)*a(2)];

% use ode45 (RK4/RK5) to solve differential equations for state X
[t,X] = ode45(f,[0 tmax],[1 1 1]); 

close all
plotLorenzSolution(X, t);

% use differential equations to obtain derivatives of state X for all t in
% the sample
for i = 1:size(t,1)
  dX(i,:) = f(t(i), X(i,:));
end

% add noise to the derivative data
eta = 0.5;      % noise magnitude
dX = dX + eta * randn(size(dX));


%% SINDy Algorithm:

% library of nonlinear functions 
lib = @(t,a) [1; a(1); a(2); a(3); a(1)^2; a(2)^2; a(3)^2; a(1)*a(2); a(1)*a(3); ...
        a(2)*a(3); a(1)^3; a(2)^3; a(3)^3; a(1)*a(2)*a(3); a(1)^2*a(2); a(1)^2*a(3); ...
        a(2)^2*a(1); a(2)^2*a(3); a(3)^2*a(1); a(3)^2*a(2); sin(a(1)); cos(a(1)); ...
        sin(a(2)); cos(a(2)); sin(a(3)); cos(a(3))]';

% input data into library
for i = 1:size(t,1)
  Theta(i,:) = lib(i,X(i,:));
end

Xi = SINDy(Theta, dX);

% now have solution Xi that describes which nonlinear terms are active in
% the dynamics of the input system

%% Recontruct data with Xi to compare with initial solution

Fk = @(i,a) Xi' * lib(i,a)';

[t,Sol] = ode45(@(t,Sol) Fk(t,Sol),t,[1 1 1]);

ax = plotLorenzSolution(Sol, t);
plotLorenzMoving(Sol, ax);
% 
