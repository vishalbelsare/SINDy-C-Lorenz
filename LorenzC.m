function f = LorenzC (t, a)
    
    
sigma=10;
beta=8/3;
rho=28;
    % Lorenz equations
    d = randn;                   % gaussian white noise
    u = 26 - a(1) + d;
%     g = u(1);
%     U = [U;d];
%     tu = [tu,t];
    f = [-sigma*a(1) + sigma*a(2) + u; rho*a(1) - a(2) - a(1)*a(3); -beta*a(3) + a(1)*a(2)];
    % u = @(t) .5 + sin(40*t);
    % g = @(u) u ^ 3;
    
end