% compute sparse regression on dX = Theta * Xi
% regression technique used: sequential least squares
% code taken directly from Supporting Info for SINDy Paper
    
function Xi = SINDy (Theta, dXdt)

    Xi = Theta \ dXdt;

    lambda = 0.1;
    for i = 1:20
      smallinds = (abs(Xi) < lambda);
      Xi(smallinds) = 0;    % set negligible terms to 0
      for ind = 1:size(dXdt,2)   % perform regression on each vector independently
        biginds = ~smallinds(:,ind);
        Xi(biginds,ind) = Theta(:,biginds) \ dXdt(:,ind);
      end
    end
    
end