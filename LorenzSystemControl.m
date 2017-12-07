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
t = linspace(0,tmax,2000);

% Lorenz equations
d = @(t) randn; % white noise
u = @(t, a) [26 - a(1) + d(t); 0; 0];
g = @(u) u(1);
f = @(t,a) [-sigma*a(1) + sigma*a(2) + g(u(t, a)); rho*a(1) - a(2) - a(1)*a(3); -beta*a(3) + a(1)*a(2)];
% u = @(t) .5 + sin(40*t);
% g = @(u) u ^ 3;

% use ode45 (RK4/RK5) to solve differential equations for state X
[t,X] = ode45(f,t,[1 1 1]); 

close all
% plotLorenzSolution(X, t);

% use differential equations to obtain derivatives of state X for all t in
% the sample
for i = 1:size(t,1)
  dXdt(i,:) = f(t(i), X(i,:));
end
% matrix of control history
for i = 1:size(t,1)
  Y(i,:) = u(t(i), X(i,:));
end


% % add noise to the derivative data
% eta = 0.5;      % noise magnitude
% dX = dX + eta * randn(size(dx));

plotLorenzSolution(X, t);

%% SINDy Algorithm to obtain Y(X):

% library of nonlinear functions 
lib = @(t,a) [1; a(1); a(2); a(3); a(1)^2; a(2)^2; a(3)^2; a(1)*a(2); a(1)*a(3); ...
        a(2)*a(3); a(1)^3; a(2)^3; a(3)^3; a(1)*a(2)*a(3); a(1)^2*a(2); a(1)^2*a(3); ...
        a(2)^2*a(1); a(2)^2*a(3); a(3)^2*a(1); a(3)^2*a(2); sin(a(1)); cos(a(1)); ...
        sin(a(2)); cos(a(2)); sin(a(3)); cos(a(3))]';
    % need more sine terms?

% input data into library
for i = 1:size(t,1)
  Theta(i,:) = lib2(i,X(i,:),[0, 0, 0]);
end

Xiu = SINDy(Theta, Y); 



% now have solution Xi that describes which nonlinear terms are active in
% the dynamics of the input system

%% New system

% redefine control
u = @(t) [50 * sin(10 * t); 0; 0];

[t,X] = ode45(f,t,[1 1 1]); 

for i = 1:size(t,1)
  Theta(i,:) = lib2(i,X(i,:),[0,0,0]);
end

YY = Theta * Xiu;

Z = YY - Y;

for i = 1:size(t,1)
  dXdt(i,:) = f(t(i), X(i,:));
end
%% SINDYc algorithm to solve for Xi s.t. dX = Theta(X, Theta(x)*Xiu) * Xi

    

for i = 1:size(t,1)
  ThetaXY(i,:) = lib2(i,X(i,:),YY(i,:));
end

Xi = SINDy(ThetaXY, dXdt);


%% Recontruct data with Xi to compare with initial solution


% yx = @(i,a)   
Fk = @(i,a) Xi' * lib2(i,a,interp1(t,YY,i))';


[t,Sol] = ode45(@(t,Sol) Fk(t,Sol),t,[1 1 1]);

ax = plotLorenzSolution(Sol, t);
% plotLorenzMoving(Sol, ax);
% 
