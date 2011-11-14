%% load data
g= temptest('30_mars_10cp_25apr.mat');
g = [g temptest('30_mars_30cp_25apr.mat')];
g = [g temptest('30_mars_50cp_25apr.mat')];
%g = [g temptest('8mars25Arp50cp.mat')];
%g = [g temptest('8mars25arp00cp.mat')];
%g = [g temptest('run1_1-mars-2011_25arp-10cp.mat')];
%g = [g temptest('run3_1-mars-2011_25arp-30cp.mat')];
%g = [g temptest('run4_1-mars-2011_25arp-30cp.mat')];

%%
i=0;
%% checkplot
i=0;
smoyx=[];
mmoyx=[];
secrx=[];
mecrx=[];
smoyy=[];
mmoyy=[];
secry=[];
mecry=[];

summoyx=[];
sumstdx=[];
summoyy=[];
sumstdy=[];




for i=1:length(g)
    figure(1);
    mtx=[g(i).moving_trap.bead_pos_in_trap.x];
    stx=[g(i).still_trap.bead_pos_in_trap.x];
    d=g(i).moving_trap.event.appr.stop;
    f=g(i).moving_trap.event.retr_start;
    plot(stx(d:f)+mtx(d:f),'g');
    hold on;
    plot(stx(d:f),'k.');
    plot(mtx(d:f),'r.');
    %axis([0 15e4 -10 10]);
    smoyx=[smoyx mean(stx(d:f))];
    mmoyx=[mmoyx mean(mtx(d:f))];
    secrx=[secrx std(stx(d:f))];
    mecrx=[mecrx std(mtx(d:f))];


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
    
    [a,b]=hist(sty(d:f)-mean(sty(d:f)),50);
    plot(b,a,'ro--');
    hold on
    [a,b]=hist(mty(d:f)-mean(mty(d:f)),50);
    plot(b,a,'bo--');
    [a,b]=hist(sty(d:f)+mty(d:f)-mean(sty(d:f)+mty(d:f)),50);
    [mu,sig,ermu,ersig] = normfit(sty(d:f)+mty(d:f)-mean(sty(d:f)+mty(d:f)));
    mu
    sig
    ermu
    ersig
    plot(b,a,'g*--');
    fu= @(x)normpdf(x,mu,sig+err)*max(a)/normpdf(0,0,sig);
    plot(b,fu(b),'g');
    
    hold off;
    
    smoyy=[smoyy mean(sty(d:f))];
    mmoyy=[mmoyy mean(mty(d:f))];
    secry=[secry std(sty(d:f))];
    mecry=[mecry std(mty(d:f))];
    
    summoyx=[summoyx mean(stx(d:f)+mtx(d:f))];
    sumstdx=[sumstdx std(stx(d:f)+mtx(d:f))];
    summoyy=[summoyy mean(sty(d:f)+mty(d:f))];
    sumstdy=[sumstdy std(sty(d:f)+mty(d:f))];


    fprintf(' %d over %d \n',i,length(g));

    hold off;

    input('next...');
end
x=[1:length(g)];
