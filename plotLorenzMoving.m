function plotLorenzMoving(X, ax)

figure()
% clf
colormap(jet);

hold off

for i = 1:length(X(:,1))
    scatter3(X(1:i,1), X(1:i,2), X(1:i,3), 20, 1:i);
    axis([ax.XLim(1) ax.XLim(2) ax.YLim(1) ax.YLim(2) ax.ZLim(1) ax.ZLim(2)]);
    view(-130,10);
    grid off
    pause(.000001);
end

end