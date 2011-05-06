%% load data
g= temptest('30_mars_10cp_25apr.mat');
g = [g temptest('30_mars_30cp_25apr.mat')];
g = [g temptest('30_mars_50cp_25apr.mat')];
g = [g temptest('8mars25Arp50cp.mat')];
g = [g temptest('8mars25arp00cp.mat')];
g = [g temptest('run1_1-mars-2011_25arp-10cp.mat')];
g = [g temptest('run3_1-mars-2011_25arp-30cp.mat')];
g = [g temptest('run4_1-mars-2011_25arp-30cp.mat')];

%%
i=0;
%% checkplot
i=0;
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


for i=1:length(g)
    figure(1);
    
    mtx=[g(i).moving_trap.bead_pos_in_trap.x];
    stx=[g(i).still_trap.bead_pos_in_trap.x];
    d=g(i).moving_trap.event.appr.stop;
    f=g(i).moving_trap.event.retr_start;
    n=floor(sqrt(f-d));
    nn= [nn n];
    plot(stx(d:f)+mtx(d:f),'g');
    hold on;
    plot(stx(d:f),'k.');
    plot(mtx(d:f),'r.');
    %axis([0 15e4 -10 10]);
    s.moy.x=[s.moy.x mean(stx(d:f))];
    m.moy.x=[m.moy.x mean(mtx(d:f))];
    s.ecr.x=[s.ecr.x std(stx(d:f))];
    m.ecr.x=[m.ecr.x std(mtx(d:f))];


    hold off;
    
    figure(2);
    hold off;
    mty=[g(i).moving_trap.bead_pos_in_trap.y];
    sty=[g(i).still_trap.bead_pos_in_trap.y];
    d=g(i).moving_trap.event.appr.stop;
    f=g(i).moving_trap.event.retr_start;
    plot(sty(d:f)+mty(d:f),'+g');
    hold on;
    plot(sty(d:f),'k.');
    plot(mty(d:f),'r.');
    hold off
    
    figure(3)
    
    [a,b]=hist(sty(d:f)-mean(sty(d:f)),n);
    plot(b,a,'r.');
    hold on
    [a,b]=hist(mty(d:f)-mean(mty(d:f)),n);
    plot(b,a,'b.');
    [a,b]=hist(sty(d:f)+mty(d:f)-mean(sty(d:f)+mty(d:f)),n);
    [mu,sig,ermu,ersig] = normfit(sty(d:f)+mty(d:f)-mean(sty(d:f)+mty(d:f)));

    plot(b,a,'g*--');
    fu= @(x)normpdf(x,mu,sig)*max(a)/normpdf(0,0,sig);
    plot(b,fu(b),'k--','LineWidth',1.2);
    err= [err sum((a-fu(b)).^2)/length(mty(d:f))];
    
    hold off;
    
    s.moy.y=[s.moy.y mean(sty(d:f))];
    m.moy.y=[m.moy.y mean(mty(d:f))];
    s.ecr.y=[s.ecr.y std(sty(d:f))];
    m.ecr.y=[m.ecr.y std(mty(d:f))];
    
    somme.moy.x=[somme.moy.x mean(stx(d:f)+mtx(d:f))];
    somme.std.x=[somme.std.x std(stx(d:f)+mtx(d:f))];
    somme.moy.y=[somme.moy.y mean(sty(d:f)+mty(d:f))];
    somme.std.y=[somme.std.y std(sty(d:f)+mty(d:f))];


    fprintf('\n%d over %d ',i,length(g));

    hold off;

    %input('next...');
end
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
