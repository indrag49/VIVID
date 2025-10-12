function impactOscillator

lw=3;

t = linspace(-3.4, 44.3, 500);
x = cos(t) + 0.4*t;
y = sin(t)-3;

figure(3)
set(gcf,'Position',[500,200,840,420])
set(gca,'Position',[0,0,1,1])
clf
hold on
plot(x, y, 'k-', 'LineWidth', lw)

%spring tails
plot([-15 -2.36], [-2.8 -2.8], 'k-', 'LineWidth', lw)
plot([18.53 31], [-2.7 -2.7], 'k-', 'LineWidth', lw)

% plot([-15, -15], [7, -10], 'k-', 'LineWidth', lw) %left wall

vertices = [
   -15, -8;
   -15,   5;
   -20,   5;
   -20,  -8;
   -15, -8
];
xv = vertices(:,1);
yv = vertices(:,2);

patch(xv, yv, [0.44, 0.50, 0.56], 'EdgeColor','none');

plot([31, 31, 40, 40, 31], [.2, -3.2, -3.2, .2, .2], 'k-', 'LineWidth', 5) %block

vertices = [
   31, .2;
   31, -3.2;
   40, -3.2;
   40, .2;
   31, .2
];
xv = vertices(:,1);
yv = vertices(:,2);

patch(xv, yv, validatecolor("#bebebe"), 'EdgeColor','none');

plot([10, 10, 16, 16, 10], [-1.4, 0, 0, -1.4, -1.4], 'k-', 'LineWidth', lw)
plot([7, 10], [-1.4, -1.4], 'k-', 'LineWidth', lw)
plot([7, 10], [0, 0], 'k-', 'LineWidth', lw)
plot([-15, 10], [-.7, -.7], 'k-', 'LineWidth', lw)
plot([16, 31], [-.7, -.7], 'k-', 'LineWidth', lw)

% plot([45, 45], [5, -8], 'k-', 'LineWidth', lw) %right wall

vertices = [
   45, 5;
   45, -8;
   50, -8;
   50, 5
];
xv = vertices(:,1);
yv = vertices(:,2);

patch(xv, yv, validatecolor("#708090"), 'EdgeColor','none');

plot([31, 40], [1., 1.], 'k-', 'LineWidth', lw)
plot([38, 40], [1.4, 1.], 'k-', 'LineWidth', lw)
plot([38, 40], [0.6, 1.], 'k-', 'LineWidth', lw)
plot([31, 33], [1., 1.4], 'k-', 'LineWidth', lw)
plot([31, 33], [1., 0.6], 'k-', 'LineWidth', lw)

text(33.5, 1.8, '$x(t)$', 'Interpreter','latex','FontSize',25)

plot([22, 27.5], [-1.5, -1.5], 'k-', 'LineWidth', 5)
plot([27.5, 27.5, 30.8, 27.5], [-1, -2, -1.5, -1], 'k-', 'LineWidth', lw)

vertices = [
   27.5, -1;
   27.5, -2;
   30.8, -1.5
];
xv = vertices(:,1);
yv = vertices(:,2);

patch(xv, yv, 'black', 'EdgeColor','none');

text(33, -4., 'block', 'Interpreter','latex','FontSize',25)

plot([26, 24.33], [-1.5, -6], 'k-', 'LineWidth', 1)
text(20, -6.4, '$\mathcal{A} \cos(\omega t)$', 'Interpreter','latex','FontSize',25)

text(37, -7.5, '$x=0$', 'Interpreter','latex','FontSize',25)
text(50.5, -1.5, 'wall', 'Interpreter','latex','FontSize',25)

text(11, 0.7, '$\zeta$', 'Interpreter','latex','FontSize',25)

plot([45, 41], [-5, -7], 'k-', 'LineWidth', 1)

axis("off")

