clear all
fs = 16; % font size
lw = 2;  % line wideth 
ms = 16; % marker size

m  = 75;
x  = linspace(-pi,pi,m);
y1 = cosh(x);
y2 = sinh(x);

figure
plot( x,y1,'bx-', 'lineWidth',lw, 'MarkerSize',ms  );
hold on
plot( x,y2,'ro--', 'lineWidth',lw, 'MarkerSize',ms  );
hold off
xlabel( 'x' );
ylabel( 'y' );
legend( 'cosh(x)','sinh(x)', 'Location','SouthEast' );
set(gca,'FontSize',fs);
plotName = sprintf('images/plotExample.eps');
fprintf('Saving file=[%s]\n',plotName);
print('-depsc2',plotName);
