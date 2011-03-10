function r = testfit()
    %%
    d=1:0.1:100;
    r=[];
    for i=1:100
        %%
        t.dr=normrnd(80,10);
        t.fr=normrnd(1,1);
        t.er=normrnd(1,1);
        yorig=younghertz(d,t.dr,t.fr,t.er)+normrnd(0,20,size(d));
        %%
        [D,F,E,err]= youngfit(d,yorig);
        %%
        t.dg=D;
        t.fg=F;
        t.eg=E;
        t.err=err; 
        r=[r,t];
        fprintf('%d/100, err,%d\n',i,t.err);
        figure(1);
        hold on 
        plot(d,yorig,'x');
        plot(d,younghertz(d,D,F,E),'r*');
        hold on
    end
end
