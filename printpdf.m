function printpdf(i)
figure(i)
set(gcf,'paperunits','centimeters')
set(gcf,'papersize',[29.7,21]) % Desired outer dimensions
set(gcf,'paperposition',[0,9,29,20]);
print -dpdf myfigure.pdf % Place plot on figure