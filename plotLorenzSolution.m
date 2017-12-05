function ax = plotLorenzSolution(X, t)


figure()
clf
scatter3(X(:,1), X(:,2), X(:,3), 6, t);
colormap(summer);
view(-130,10);
ax = gca;
hold off

end