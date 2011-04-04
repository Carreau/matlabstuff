function r = testfit()
    %%
    % dinit=8;
    %Einit=3e-12;
    %f0init=0;
    %Drift = 0 ;
    d=1:0.1:45;
    r=[];
    fprintf('\niteration:       \n\n\r');
    nmax=5000;
    for i=1:nmax
        %%
        t.dr=normrnd(15,6);       %touch distance
        t.fr=normrnd(3e-12,2E-12); 
        t.er=normrnd(4e-12,3e-12); % young modulus
        t.drift = normrnd(-3E-12,3E-12);
        yorig=younghertz(d,t.dr,t.fr,t.er,t.drift)+1*normrnd(0,20,size(d))*3e-12;
        %%
        [D,F,E,err,drift]= youngfit(d,yorig);
        %%
        t.dg=D;
        t.fg=F;
        t.eg=E;
        t.errg=err;
        t.driftg=drift;
        r=[r,t];
        fprintf('\b\b\b\b\b\b\b\b\b');
        fprintf('%04d/%04d',i,nmax);
        %figure(1);
        hold on 
        %plot(d,yorig,'x');
        %plot(d,younghertz(d,D,F,E,drift),'r*');
        hold on
    end
    %figure(2);
    %clf
    %plot([r.er],[r.eg],'+');
    %title('young modulus');
    
    %figure(3);
    %plot([r.dr],[r.dg],'+');
    %title('touch distance');
    
    %figure(4);
    %dd =[[r.dr]-[r.dg],[r.er]-[r.eg]];
    %dd = reshape(dd,[length(dd)/2,2]);
    %hist3(dd,[20,20]);
    %set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
    c = clock;
    str = sprintf('stat_%i_%i-%i-%i_%ih%i.mat',nmax,c(1),c(2),c(3),c(4),c(5));
    save(str,'r');
end
