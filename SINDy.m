% compute sparse regression on dX = Theta * Xi
% regression technique used: sequential least squares
% code taken directly from Supporting Info for SINDy Paper
    
function Xi = SINDy (Theta, dX)

    Xi = Theta \ dX;

    lambda = 0.1;
    for i = 1:20
      smallinds = (abs(Xi) < lambda);
      Xi(smallinds) = 0;    % set negligible terms to 0
      for ind = 1:size(dX,2)   % perform regression on each vector independently
        biginds = ~smallinds(:,ind);
        Xi(biginds,ind) = Theta(:,biginds) \ dX(:,ind);
      end
    end
    disp(size(Xi));
    
end