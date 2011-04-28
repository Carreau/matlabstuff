function [d0,E,err,flag] = yf(d,yorig,Starting)
    %%
    %%
    [Estimates,~,flag]=fminsearch(@ye,Starting,[],d,yorig);
    %Estimates = Starting;
    d0 = Estimates(1);
    E  = Estimates(2);
    err=ye(Estimates,d,yorig);
    if(flag)
        figure(1)
        %hold off;
        plot(d,yorig,'+');
        hold on;
        plot(d,yz(d,d0,E),'r*');
    end
    %hold off;
    %printpdf(1,str);
end
     

