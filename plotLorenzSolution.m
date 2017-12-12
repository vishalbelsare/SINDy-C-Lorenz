function ax = plotLorenzSolution(X, t)
x = X(:,1);
y = X(:,2);
z = X(:,3);

figure()
clf
% C = (x + y - z) / (x + y + z);
plot3(X(:,1), X(:,2), X(:,3));
% colormap(jet);
fontsize = 20;
xlabel('$x$','Interpreter','Latex','Fontsize',fontsize)
ylabel('$y$','Interpreter','Latex','Fontsize',fontsize)
zlabel('$z$','Interpreter','Latex','Fontsize',fontsize)
view(-130,10);
ax = gca;
hold off

end