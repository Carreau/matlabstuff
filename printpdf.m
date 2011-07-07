function printpdf(i,filename)
    % use : printpdf(figurenumber, filename')
    figure(i)
    set(gcf,'paperunits','centimeters')
    set(gcf,'papersize',[29.7,21]) % Desired outer dimensions
    set(gcf,'paperposition',[0,0,29,20]);
    
    strstr = sprintf('%s-%s',datestr(now,'YYYY-mmm-dd'),filename);
    print('-dpdf',strstr); % Place plot on figure
    %plot2svg([297,110],i);

end
