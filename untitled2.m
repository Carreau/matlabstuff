%% load data
g= temptest('30_mars_10cp_25apr.mat');
g = [g temptest('30_mars_30cp_25apr.mat')];
g = [g temptest('30_mars_50cp_25apr.mat')];
g = [g temptest('8mars25Arp50cp.mat')];
%g = [g temptest('8mars25arp00cp.mat')];
%g = [g temptest('run1_1-mars-2011_25arp-10cp.mat')];
%g = [g temptest('run3_1-mars-2011_25arp-30cp.mat')];
%g = [g temptest('run4_1-mars-2011_25arp-30cp.mat')];

%%
i=0;

%% check if sum of force is constent during the resting
clear i s m somme
s.moy.x=[];
m.moy.x=[];
s.ecr.x=[];
m.ecr.x=[];
s.moy.y=[];
m.moy.y=[];
s.ecr.y=[];
m.ecr.y=[];

somme.moy.x=[];
somme.std.x=[];
somme.moy.y=[];
somme.std.y=[];
err=[];
nn=[];
clear f;
tic;
parfor i=1:length(g)
    figure(1);
    
    
    mtx=[g(i).moving_trap.bead_pos_in_trap.x];
    stx=[g(i).still_trap.bead_pos_in_trap.x];
    d=g(i).moving_trap.event.appr.stop;
    f=g(i).moving_trap.event.retr_start;
    n=floor(sqrt(f-d));
    range=[d:f];
    nn= [nn n];
    %plot(stx(d:f)+mtx(d:f),'g');
    %hold on;
    %plot(stx(d:f),'k.');
    %plot(mtx(d:f),'r.');
    %axis([0 15e4 -10 10]);
    s(i).moy.x=mean(stx(range));
    m(i).moy.x=mean(mtx(range));
    s(i).ecr.x=std(stx(range));
    m(i).ecr.x=std(mtx(range));

    
    hold off;
    
    figure(2);
    hold off;
    mty=[g(i).moving_trap.bead_pos_in_trap.y];
    sty=[g(i).still_trap.bead_pos_in_trap.y];
    d=g(i).moving_trap.event.appr.stop;
    f=g(i).moving_trap.event.retr_start;
    %plot(sty(d:f)+mty(d:f),'+g');
    hold on;
    %plot(sty(d:f),'k.');
    %plot(mty(d:f),'r.');
    hold off
    
    figure(3)
    
    [a,b]=hist(sty(d:f)-mean(sty(d:f)),n);
    %plot(b,a,'r.');
    %hold on
    [a,b]=hist(mty(d:f)-mean(mty(d:f)),n);
    %plot(b,a,'b.');
    [a,b]=hist(sty(d:f)+mty(d:f)-mean(sty(d:f)+mty(d:f)),n);
    [mu,sig,ermu,ersig] = normfit(sty(d:f)+mty(d:f)-mean(sty(d:f)+mty(d:f)));

    %plot(b,a,'g*--');
    fu= @(x)normpdf(x,mu,sig)*max(a)/normpdf(0,0,sig);
    %plot(b,fu(b),'k--','LineWidth',1.2);
    err= [err sum((a-fu(b)).^2)/length(mty(d:f))];
    
    hold off;
    
    s(i).moy.y=mean(sty(d:f));
    m(i).moy.y=mean(mty(d:f));
    s(i).ecr.y=std(sty(d:f));
    m(i).ecr.y=std(mty(d:f));
    
    somme(i).moy.x=mean(stx(d:f)+mtx(d:f));
    somme(i).std.x=std(stx(d:f)+mtx(d:f));
    somme(i).moy.y=mean(sty(d:f)+mty(d:f));
    somme(i).std.y=std(sty(d:f)+mty(d:f));
    

    fprintf('\n%d over %d ',i,length(g));

    hold off;

    %input('next...');
end
clear n nn range 
toc
x=[1:length(g)];




%% %%
cps =[g(garde==1).applyfun(@(x)x.cp)]; taus= [g(garde==1).applyfun(@(x)x.fitvalue.estimates.tau)];
figure(4);
want = (abs(taus) < 5);
cps=cps(want);
taus=taus(want);
hold off;
clf;
%plot(cps,taus,'+');
hold on ;
eerrw=max(err(want))-err(want);
%eerrw=1./err(want);

for vcp=[0 10 30 50]
    %select datapoint by cp
    dp = taus(cps == vcp);
    errw = eerrw(cps == vcp);
    %errorbar(vcp+2,mean(dp),std(dp),'ro');
    %
    tmp_med=median(dp);
    tmp_std=5*std(dp);
    strip = (dp > tmp_med-tmp_std) & (dp < tmp_med+tmp_std );
    sum(strip)/length(errw)
    dp=dp(strip);
    errw=errw(strip);
    plot(vcp*ones(size(dp)),dp,'+');
    %errorbar(vcp-1,tmp_med,tmp_std,'ro');
    
    errorbar(vcp+1,wmean(dp,errw),sqrt(var(dp,errw)),'ro');
    axis([-1 55 -0.03 0.25]);
end

%%
garde=[];
for i=1:length(g)
   err(i);
   g(i).showfitTau;
   inp=input('garder ? [y/N] ','s');
   if(inp=='n')
       garde=[garde 0];
       disp('on garde pas');
   else
       garde=[garde 1];
       disp('on garde');
   end
end

%% 
mf=[];
figure(1);
hold on;
imax = length(f)
for i=1:imax
    exp=f(i);
    start=1;
    stop = exp.moving_trap.event.appr.stop;
    %d= exp.bead_distance(start:stop);
    force= exp.still_trap_force.tangent(start:stop);
    %plot(d,f);
    fprintf('%d/%d\n',i,imax);
    %input('next...');
    mf = [mf max(force)];
end

%%
dd = [];
for i=1:length(f)
    exp=f(i);
    start=1;
    stop = exp.moving_trap.event.appr.stop;
    d= exp.bead_distance(start:stop);
    force= exp.still_trap_force.tangent(start:stop);
    ee=abs(force-1e-11);
    dd=[dd d(find( ee == min(ee)))];
end

%%
cps = [f.cp];
figure(1)
hold on;
for vcp=[0 10 30 50]
    dp = dd(cps == vcp);
    plot(vcp*ones(size(dp)),dp,'+');
    errorbar(vcp+1,mean(dp),std(dp),'ro');
end

%% 
powercoeff=[];
figure(1);
clf; hold on;
figure(2)
clf;
for i=1:length(f)
    Ehair=[];
    Dhair=[];
    exp=f(i);
    start=1;
    stop = exp.moving_trap.event.appr.stop;
    d= exp.bead_distance(start:stop);
    force= abs(exp.still_trap_force.tangent(start:stop));
    n=50;m=10;
    schunk=floor((stop-start)/n);
    lchunk=m*schunk;
    figure(3);
    clf;
    hold on;
    plot(d,force,'+');
    ppl=[];
    for j=1:(n-m)
         tempdd = mean(d(j*schunk:j*schunk+lchunk));
         Dhair = [Dhair tempdd];
         deltad = d(j*schunk)-d(j*schunk+lchunk);
         deltaf = force(j*schunk)-force(j*schunk+lchunk);
         p=polyfit(d(j*schunk:j*schunk+lchunk),force(j*schunk:j*schunk+lchunk),1);
         ppl= [ppl abs(p(1)*tempdd)];
         clear p;
         plot([d(j*schunk) d(j*schunk+lchunk)],[force(j*schunk) force(j*schunk+lchunk)],'ro-');
         Ehair = [Ehair abs(deltaf*tempdd/deltad)];
         figure(3)
         clear tempdd
    end
    figure(2);
    clf;
    %semilogy(Dhair,Ehair,'+-');
    plot(log(Dhair),log(Ehair),'+-');
    hold on;
    %semilogy(Dhair,ppl,'ro-');
    plot(log(Dhair),log(ppl),'ro-');
    l = length(Dhair);
    nbr=10;
    p = polyfit(log(Dhair(l-20:l)),log(ppl(l-20:l)),1);
    powercoeff = [powercoeff p(1)];
    plot(log(Dhair),p(1)*log(Dhair)+p(2),'g--');
    plot(log(Dhair),-2*log(Dhair)-20,'k--');
    %isd = 1 ./ (Dhair .^2)
    %semilogy(Dhair,isd,'.k');
    figure(1);
    clf;
    plot(powercoeff,'+');
    %clear p;
    fprintf('%d/%d',i,length(f));
    %input('next...');
end
