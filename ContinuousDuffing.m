clear


m = 10000;
tmax = 80;
t = linspace(0,tmax,m);

alpha = 1;
beta = 5;
gamma = 8;% d = 0.1 not in chaotic regime, d = bif variable
delta = .2;
omega = 1.4; 

% a2. = a1;
u = @(t) gamma*cos(omega*t);
f = @(t,a) [a(2); -(alpha*a(1) + beta*a(1)^3 + delta*a(2)) + u(t)];

[t,X] = ode45(f, t, [1 0]);

figure(1)
clf
plot(X(:,1), X(:,2),'k');
saveas(gcf,'CD_A1.png');

% compute derivatives of the state X for all t in the sample
% and measure the value of the forcing at each sample
for i = 1:m
  dXdt(i,:) = f(t(i), X(i,:));
  Y(i) = u(t(i));
end

% add noise to the derivative data
eta = 0.1;      % noise magnitude
dXdt = dXdt + eta * randn(size(dXdt));




%% SINDy Algorithm:

% library of nonlinear functions including control
lib = @(t,a,b) [1; a(1); a(2); b; a(1)^2; a(2)^2; a(1)*a(2); a(1)*b; a(2)*b; ...
        a(1)^3; a(2)^3; a(1)^2*a(2); a(2)^2*a(1); a(1)^2*b; a(1)*b^2; a(2)^2*b; a(2)*b^2; ...
        sin(a(1)); cos(a(1)); sin(a(2)); cos(a(2)); sin(b); cos(b)]';
    

% input data into library
for i = 1:m
  Theta(i,:) = lib(i,X(i,:),Y(i));
end

Xi = SINDy(Theta, dXdt);

% now have solution Xi that describes which nonlinear terms are active in
% the dynamics of the input system


%% Recover data with Xi 


Fk = @(i,a,b) Xi' * lib(i,a,b)';

[t,Xs] = ode45(@(t,Xs) Fk(t,Xs,u(t)), t, [1 0]);
figure(2)
clf
plot(Xs(:,1), Xs(:,2),'k');
saveas(gcf,'CD_S1.png');


gamma = 10;
u = @(t) gamma*cos(omega*t);
f = @(t,a) [a(2); -(alpha*a(1) + beta*a(1)^3 + delta*a(2)) + u(t)];

[t,Xs1] = ode45(@(t,Xs1) Fk(t,Xs1,u(t)), t, [1 0]);
[t,X1] = ode45(f, t, [1 0]);


figure(3)
clf
plot(Xs1(:,1), Xs1(:,2),'k');
saveas(gcf,'CD_A2.png');
figure(4)
clf
plot(X1(:,1), X1(:,2),'k');
saveas(gcf,'CD_S2.png');

% Error = SimulationError(X(:,2), Xs(:,2));


