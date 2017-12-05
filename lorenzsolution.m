
clear
% Solve the Lorenz system in the time interval [0,100] with initial conditions [1,1,1]
% to obtain the state data X
%
sigma=10;
beta=8/3;
rho=28;
tmax = 10;
%
f = @(t,a) [-sigma*a(1) + sigma*a(2); rho*a(1) - a(2) - a(1)*a(3); -beta*a(3) + a(1)*a(2)];
%'f' is the set of differential equations and 'a' is an array containing values of x,y, and z variables.
%'t' is the time variable
%
[t,X] = ode45(f,[0 tmax],[1 1 1]); % 'ode45' uses adaptive Runge-Kutta method of 4th and 5th order to solve differential equations
%
% disp(X(1,:));
close all

% disp(size(X));

% ax = plotLorenzSolution(X, t);

% plotLorenzMoving(X, ax);


% build time derivatives based on Lorenz equations
for i = 1:size(t,1)
  dX(i,:) = f(t(i), X(i,:));
end
%
% add noise to the derivatives
eta = .5; % noise magnitude
dX = dX + eta*randn(size(dX));
%
% build library of possible functions
libT = @(t,a) [1; a(1); a(2); a(3); a(1)^2; a(2)^2; a(3)^2; a(1)*a(2); a(1)*a(3); ...
        a(2)*a(3); a(1)^3; a(2)^3; a(3)^3; a(1)*a(2)*a(3); a(1)^2*a(2); a(1)^2*a(3); ...
        a(2)^2*a(1); a(2)^2*a(3); a(3)^2*a(1); a(3)^2*a(2); sin(a(1)); cos(a(1)); ...
        sin(a(2)); cos(a(2)); sin(a(3)); cos(a(3))]';
%
for i = 1:size(t,1)
  Theta(i,:) = libT(i,X(i,:));
end
%
%
% compute sparse regression: sequential least squares, from Appendix B
Xi = Theta \ dX;
%
lambda = 0.1;
for i = 1:20
  smallinds = (abs(Xi) < lambda);
  Xi(smallinds) = 0;    % set negligible terms to 0
  for ind = 1:size(X,2)   % perform regression on each vector independently
    biginds = ~smallinds(:,ind);
    Xi(biginds,ind) = Theta(:,biginds) \ dX(:,ind);
  end
end

% need to sum column vectors & return in column vector
% for i = 1:size(t,1)
%   Theta(i,:) = libT(i,X(i,:));
% end



% function F = Fk(t, i, a) 
% FK = @(i,a) NewTheta = libT(i,a(:));

Fk = @(i,a) Xi' * libT(i,a)';


[t,Sol] = ode45(@(t,Sol) Fk(t,Sol), t,[1 1 1]);

disp(size(Sol));

ax = plotLorenzSolution(Sol, t);
plotLorenzMoving(Sol, ax);
% 
























