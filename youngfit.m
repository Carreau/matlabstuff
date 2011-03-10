function [d0,f0,E,err] = youngfit(d,yorig)
    %%
    dinit=12;
    Einit=3e-12;
    f0init=0;
    Starting = [dinit,f0init,Einit];
    %%
    Estimates=fminsearch(@youngerr,Starting,[],d,yorig);
    d0 = Estimates(1);
    f0 = Estimates(2);
    E  = Estimates(3);
    %%
    err=youngerr(Estimates,d,yorig);
    figure(1)
    hold off;
    plot(d,yorig,'x');
    hold on;
    plot(d,younghertz(d,d0,f0,E),'r*');
    hold off;
end
     

