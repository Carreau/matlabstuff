function printpdf(i,filename)
    figure(i)
    set(gcf,'paperunits','centimeters')
    set(gcf,'papersize',[29.7,21]) % Desired outer dimensions
    set(gcf,'paperposition',[0,0,29,20]);
    print('-dpdf',filename) % Place plot on figure
end