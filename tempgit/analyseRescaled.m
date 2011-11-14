%%% analysed rescales curves.
%% laod the curves in "q" %%

%q= load('rescalesvalues.mat');
q= load('filterd_master12.mat');


%% Transform it from structure of array to array of structure
% (simpler for filtering)

len = length(q.dmax_a);
r=zeros(0,len);
for i=1:len
    r(i).dmax   = q.dmax_a(i);
    r(i).ddemi  = q.ddemi_a(i);
    r(i).dquart = q.dquart_a(i);
    r(i).cp     = q.cp_a(i);
    r(i).fmax   = q.fmax_a(i);
end 
clear j;
    
%% first remove all the curves which have a too small max force.
w= [r(:).fmax] < 0.5e-11;
r(w)=[];
fprintf('there is %d experiment left after force filtering\n',length(r));

%% filtre par demidistance >15um
w= [r(:).ddemi] > 12;
r(w)=[];
fprintf('there is %d experiment left after demi distance filtering\n',length(r));

%% retrait particulier
r(11)=[];

%% now lets analyse

%% fmax
figure(1)
plot([r.fmax],'+')
title('liste of max forces')

figure(2)
plot([r.cp],[r.fmax],'r+')
title('liste of max forces, fonction of caping concentration');
xlabel('concentration in caping protéine');
ylabel('max force');
xl=xlim; xlim([xl(1)-1 xl(2)+1]);

%% dmax
ofset=0.0
figure(3)
plot(([r.dmax]-ofset),'+')
title('liste of dmax')

figure(4)
plot([r.cp],([r.dmax]-ofset),'r+')
title('liste of max distances, fonction of caping concentration');
xlabel('concentration in caping protéine');
ylabel('max distance');
xl=xlim; xlim([xl(1)-1 xl(2)+1]);

%% ddemi

figure(11)
plot(([r.dquart]-ofset),'+')
title('liste of demi distance')

figure(12)
plot([r.cp],([r.dquart]-ofset),'r+')
title('liste of "quart distance", fonction of caping concentration');
xlabel('concentration in caping protéine');
ylabel('quart distance');
xl=xlim; xlim([xl(1)-1 xl(2)+1]);

%% dquart

figure(5)
plot(([r.ddemi]-ofset),'+')
title('liste of demi distance')

figure(6)
plot([r.cp],([r.ddemi]-ofset),'r+')
title('liste of "demi distance", fonction of caping concentration');
xlabel('concentration in caping protéine');
ylabel('demi distance');
xl=xlim; xlim([xl(1)-1 xl(2)+1]);

%% fmax/ddemi
figure(7)
plot([r.fmax] ./ ([r.ddemi]-ofset),'+')
title('liste of fmax/ddemi')

figure(8)
plot([r.cp],[r.fmax] ./ ([r.ddemi]-ofset),'r+')
title('liste of "fmax/ddemi", fonction of caping concentration');
xlabel('concentration in caping protéine');
ylabel('fmax/ddemi');
xl=xlim; xlim([xl(1)-1 xl(2)+1]);


%% fmax*ddemi
figure(9)
plot([r.fmax] .* ([r.ddemi] -ofset),'+')
title('fmax*ddemi')

figure(10)
plot([r.cp],[r.fmax] .* ([r.ddemi] -ofset),'r+')
title('"fmax*ddemi", vs [cp]');
xlabel('concentration in caping protéine');
ylabel('fmax*ddemi');
xl=xlim; xlim([xl(1)-1 xl(2)+1]);

%% 
figure(31)
clf
figure(30)
clf

% essaie sur la variance 
for ccp=[0 10 30 50];
    w = [r.cp] == ccp;
    s=r(w);
    vvar = @(x,y) var([y.fmax] .* ([y.ddemi]-x));%/mean([y.fmax] ./ ([y.ddemi]-x));
    svvar = @(x) vvar(x,s);
    asvvar = @(x) arrayfun(svvar,x);
    interval= [0.5:0.01:10.0];
    figure(30)
    hold on 
    plot(interval,asvvar(interval),'r')
    
    [minv,mini] = min(asvvar(interval));
    dmin = interval(mini);
    fprintf(' pour cp = %d variance minimale (%d) pour dofset=%d\n',ccp,minv,dmin);

    figure(31)
    hold on 
    plot(ccp,dmin,'ro');
    %for i=interval
    %   plot(i*ones(1,length(s)),[s.fmax] .* ([s.ddemi]-i),'+');
    %end
end
xl=xlim; xlim([xl(1)-1 xl(2)+1]);
yl=ylim; ylim([yl(1)-1 yl(2)+1]);

