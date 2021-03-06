clear

lib2 = @(t,a, b) [1; a(1); a(2); a(3); b(1); b(2); b(3); ...      % 1, X, Y
        a(1)^2; a(2)^2; a(3)^2; a(1)*a(2); a(1)*a(3); a(2)*a(3); ...  % X^2
        b(1)^2; b(2)^2; b(3)^2; b(1)*b(2); b(1)*b(3); b(2)*b(3); ...  % Y^2
        a(1)*b(1); a(1)*b(2); a(1)*b(3); a(2)*b(1); a(2)*b(2); ...
        a(2)*b(3); a(3)*b(1); a(3)*b(2); a(3)*b(3); ...               % X*Y
%         a(1)^3; a(2)^3; a(3)^3; a(1)*a(2)*a(3); a(1)^2*a(2); a(1)^2*a(3); ...
%         a(2)^2*a(1); a(2)^2*a(3); a(3)^2*a(1); a(3)^2*a(2);           % X^3
%         b(1)^3; b(2)^3; b(3)^3; b(1)*b(2)*b(3); b(1)^2*b(2); b(1)^2*b(3); ...
%         b(2)^2*b(1); b(2)^2*b(3); b(3)^2*b(1); b(3)^2*b(2);           % Y^3 
%         a(1)^2*b(1); a(1)^2*b(2); a(1)^2*b(3); ...
%         a(2)^2*b(1); a(2)^2*b(2); a(2)^2*b(3); ...
%         a(3)^2*b(1); a(3)^2*b(2); a(3)^2*b(3); ...      * X^2*Y
%         b(1)^2*a(1); b(1)^2*a(2); b(1)^2*a(3); ...
%         b(2)^2*a(1); b(2)^2*a(2); b(2)^2*a(3); ...
%         b(3)^2*a(1); b(3)^2*a(2); b(3)^2*a(3); ...      % Y^2*X
        sin(a(1)); cos(a(1)); sin(a(2)); cos(a(2)); sin(a(3)); cos(a(3))    % sin(X), cos(X)
        sin(b(1)); cos(b(1)); sin(b(2)); cos(b(2)); sin(b(3)); cos(b(3))    % sin(Y), cos(Y)
        ]';
%         sin(a(1)); cos(a(1)*b(1)); sin(a(2)); cos(a(2)); sin(a(3)); cos(a(3))    % sin(X), cos(X)
    


%% Set up and solve Lorenz system to get state X and dX/dt

sigma=10;
beta=8/3;
rho=28;
tmax = 20;
m = 2000;
tspan = linspace(0,tmax,m);
% global U
% global tu

% Lorenz equations
d = @(t) randn; % white noise
u = @(t, a) [26 - a(1) + 5*sin(t) + d(t); 0; 0];
g = @(U) U(1);
f = @(t,a) [-sigma*a(1) + sigma*a(2) + g(u(t, a)); rho*a(1) - a(2) - a(1)*a(3); -beta*a(3) + a(1)*a(2)];
% u = @(t) .5 + sin(40*t);
% g = @(u) u ^ 3;
% f = @(t,X) dynamics(t,X,d);

% use ode45 (RK4/RK5) to solve differential equations for state X
[tspan,X] = ode45(f, tspan, [1 1 1]); 

close all
% plotLorenzSolution(X, t);

% use differential equations to obtain derivatives of state X for all t in
% the sample
for i = 1:m
  dXdt(i,:) = f(tspan(i), X(i,:));
end
% matrix of control history
for i = 1:m
  Y(i,:) = u(tspan(i), X(i,:));
end


% % add noise to the derivative data
% eta = 0.5;      % noise magnitude
% dX = dX + eta * randn(size(dx));

plotLorenzSolution(X, tspan);

%% SINDy Algorithm to obtain Y(X):

% library of nonlinear functions 
lib = @(t,a) [1; a(1); a(2); a(3); a(1)^2; a(2)^2; a(3)^2; a(1)*a(2); a(1)*a(3); ...
        a(2)*a(3); a(1)^3; a(2)^3; a(3)^3; a(1)*a(2)*a(3); a(1)^2*a(2); a(1)^2*a(3); ...
        a(2)^2*a(1); a(2)^2*a(3); a(3)^2*a(1); a(3)^2*a(2); sin(a(1)); cos(a(1)); ...
        sin(a(2)); cos(a(2)); sin(a(3)); cos(a(3))]';
    % need more sine terms?

% input data into library
for i = 1:m
  Theta(i,:) = lib2(i,X(i,:),[0, 0, 0]);
end

Xiu = SINDy(Theta, Y); 



% now have solution Xi that describes which nonlinear terms are active in
% the dynamics of the input system

%% SINDYc algorithm to solve for Xi s.t. dX = Theta(X, Theta(x)*Xiu) * Xi

YY = Theta * Xiu;  % ok ?

for i = 1:m
  ThetaXY(i,:) = lib2(i,X(i,:),YY(i,:));
end

Xi = SINDy(ThetaXY, dXdt);


%% Recontruct data with Xi and new control to compare with initial solution


% redefine control
u = @(t) [50 * sin(10 * t); 0; 0];

% [tspan,X] = ode45(f,tspan,[1 1 1]); 

for i = 1:m
  Ynew(i,:) = u(tspan(i));
end
% 
% for i = 1:length(tspan)
%   Theta(i,:) = lib2(i,X(i,:),Ynew);
% end
% 
% Z = YY - Y;
% 
% for i = 1:length(tspan)
%   dXdt(i,:) = f(tspan(i), X(i,:));
% end
% 

Fk = @(t,a) Xi' * lib2(t,a,interp1(tspan,Ynew,t))';


[tspan,Sol] = ode45(@(t,Sol) Fk(t,Sol),tspan,[1 1 1]);

ax = plotLorenzSolution(Sol, tspan);
% plotLorenzMoving(Sol, ax);
% 
