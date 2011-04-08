function [d0,f0,E,err,Drift] = youngfit(d,yorig)
    %%
    dinit=8;
    Einit=3e-12;
    f0init=0;
    Drift = 0 ;
    Starting = [dinit,f0init,Einit,Drift];
    %%
    d=d(100:end);
    yorig=yorig(100:end);
    l=length(d);
    n=50;
    %d0      = d(1:floor(l/n):end);
    %yorig0  = yorig(1:floor(l/n):end);
    %Starting = fminsearch(@youngerr,Starting,[],d0,yorig0);
    Estimates=fminsearch(@youngerr,Starting,[],d,yorig);
    %Estimates = Starting;
    d0 = Estimates(1);
    f0 = Estimates(2);
    E  = Estimates(3);
    Drift = Estimates(4);
    err=youngerr(Estimates,d,yorig);
    
    figure(1)
    %hold off;
    plot(d,yorig,'+');
    hold on;
    plot(d,younghertz(d,d0,f0,E,Drift),'r*');
    %hold off;
    xlabel('bead distance');
    ylabel('force');
    dv = clock();
    str = sprintf('fig_%i_%i_%i-%i_%i_%i.pdf',dv(1),dv(2),dv(3),dv(4),dv(5),floor(dv(6)));
    %printpdf(1,str);
end
     

