%% load data
clear
g =    temptest('30_mars_10cp_25apr.mat');
%g = [g temptest('30_mars_30cp_25apr.mat')];
%g = [g temptest('30_mars_50cp_25apr.mat')];
%g = [g temptest('8mars25Arp50cp.mat')];
g = [g temptest('8mars25arp00cp.mat')];
g = [g temptest('run1_1-mars-2011_25arp-10cp.mat')];
g = [g temptest('run3_1-mars-2011_25arp-30cp.mat')];
g = [g temptest('run4_1-mars-2011_25arp-30cp.mat')];

%%

%% fit linéaire des courbes loglog
i=4;
f=g(i).still_trap.force.r;
d=g(i).bead_distance;
step=150;
y=f(1:step:e);
x=d(1:step:e);
step=2;
xp=x(end/4:end);
yp=y(end/4:end);

xpl=log(xp);
ypl=log(yp);

p=polyfit(xpl,ypl,1);
clf 
loglog(x(end/4:2*step:end),y(end/4:2*step:end),'b+')
hold on;
loglog(x(1:step:end/4),y(1:step:end/4),'r+');
loglog(xp(1:step:end),exp(polyval(p,xpl(1:step:end))),'-g');
q2=-21;
%p(2)
q=[-2 q2];
r=[-1.5 -22];
loglog(xp(1:step:end),exp(polyval(q,xpl(1:step:end))),'--k')
loglog(xp(1:step:end),exp(polyval(r,xpl(1:step:end))),'--y')
%ylim([0.0004    1.0000]*1e-10)
xlim([6 39])
xlabel('log(distance)')
ylabel('log(Force)')
%%
i=0;

%% check if sum of force is constent during the resting
clear i s m somme
s.moy.x=[];
m.moy.x=[];
s.ecr.x=[];
m.ecr.x=[];
s.moy.y=[]
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
	%figure(1);
	
	mtx=[g(i).moving_trap.bead_pos_in_trap.x];
	stx=[g(i).still_trap.bead_pos_in_trap.x];
	d=g(i).moving_trap.event.appr.stop;
	f=g(i).moving_trap.event.retr_start;
	mtx=mtx(d:f);
	stx=stx(d:f);
	n=floor(sqrt(f-d));
	nn= [nn n];
	%plot(stx+mtx,'g');
	%hold on;
	%plot(stx,'k.');
	%plot(mtx,'r.');
	%axis([0 15e4 -10 10]);
	s(i).moy.x=mean(stx);
	m(i).moy.x=mean(mtx);
	s(i).ecr.x=std(stx);
	m(i).ecr.x=std(mtx);

	
	%hold off;
	
	%figure(2);
	%hold off;
	mty=[g(i).moving_trap.bead_pos_in_trap.y];
	sty=[g(i).still_trap.bead_pos_in_trap.y];
	mty=mty(d:f);
	sty=sty(d:f);
    
	%plot(sty+mty,'+g');
	hold on;
	%plot(sty,'k.');
	%plot(mty,'r.');
	hold off
	
	%figure(3)
	
	[a,b]=hist(sty-mean(sty),n);
	%plot(b,a,'r.');
	%hold on
	[a,b]=hist(mty-mean(mty),n);
	%plot(b,a,'b.');
	[a,b]=hist(sty+mty-mean(sty+mty),n);
	[mu,sig,ermu,ersig] = normfit(sty+mty-mean(sty+mty));

	%plot(b,a,'g*--');
	fu= @(x)normpdf(x,mu,sig)*max(a)/normpdf(0,0,sig);
	%plot(b,fu(b),'k--','LineWidth',1.2);
	err= [err sum((a-fu(b)).^2)/length(mty)];
	
	hold off;
	
	s(i).moy.y=mean(sty);
	m(i).moy.y=mean(mty);
	s(i).ecr.y=std(sty);
	m(i).ecr.y=std(mty);
	
	somme(i).moy.x=mean(stx+mtx);
	somme(i).std.x=std(stx+mtx);
	somme(i).moy.y=mean(sty+mty);
	somme(i).std.y=std(sty+mty);
	

	fprintf('\n%d over %d ',i,length(g));

	hold off;

	%input('next...');
end
clear n nn range 
toc
x=[1:length(g)];




%% %% tracer de donner à revoir
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

%% demande si garder ou pas la courbe (a la main)
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
%% f==g 
f=g;

%% liste des forces maximales
mf=[];
%figure(1);
%hold on;
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
clear exp start stop imax force i

%%
f=g(mf > 0.8e-11);

%% extract 'touch point' @tf (N)
tf = 0.8e-11;
imax = length(f);
dd = zeros(size(imax));
parfor i=1:imax
	exp=f(i);
	start=1;
	stop = exp.moving_trap.event.appr.stop;
	d= exp.bead_distance(start:stop);
	force= exp.still_trap_force.tangent(start:stop);
	ee=abs(force-tf);
	dd(i)=d( ee == min(ee));
    fprintf('%d/%d\n',i,imax);
end
clear exp start stop d force ee imax;

%% plot 
cps = [f.cp];
tm = [f.time_m];
figure(1)
clf;
hold on;
for vcp=[0 10 30 50]
	dp = dd(cps == vcp);
	ttm = tm(cps == vcp);
	xx = (1:(length(ttm)))-length(ttm)/2;
	ttm2 = [xx;ttm];
	ttm2 = sortrows(ttm2',2)';
	plot(vcp*ones(size(dp))+(ttm2(2,:)-(min(ttm2(2,:))+max(ttm2(2,:)))/2)/3,dp,'+');
	errorbar(vcp,mean(dp),std(dp),'ro');
end
xlabel('caping + ofset time');
ylabel('touch point');
title('touch point as a fonction of caping concentrationand time');

%% try to fit power coefficient 
powercoeff=[];errcoeff=[];
powercoeff2=[];
figure(1);
clf; hold on;
figure(2)
clf;
i=1;
ldh0=0;
leh0=0;
for i=1
    
	Ehair=[];
	Dhair=[];
	exp=f
	start=1;
	stop = exp.moving_trap.event.appr.stop;
	d= exp.bead_distance(start:stop);
	force= abs(exp.still_trap_force.tangent(start:stop));
    %d=[10:-0.001:1.3];
    %fc=@(x) 1/(x-i/10)^2;
    %force = arrayfun(fc,d)+randn(size(d))/10;
    %d= d+randn(size(d))/10;
    start = 1;
    stop = length(d);
    %
	n=450;m=50;
	schunk=floor((stop-start)/n);
	lchunk=m*schunk;
	figure(3);
	clf;
	hold on;
	plot(d,force,'+');
	ppl=[];
    ldh0=0;
    
	for j=10:(n-m)
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
    if(i==1)
        ldh0=0;log(Dhair(end));
        leh0=0;log(Ehair(end));
        %continue
    end
	figure(2);
	clf;
	semilogy(Dhair,Ehair,'+-');
	plot(log(Dhair)-ldh0,log(Ehair)-leh0,'+-');
	hold on;
	%semilogy(Dhair,ppl,'ro-');
	plot(log(Dhair)-ldh0,log(ppl)-leh0,'ro-');
	l = length(Dhair);
	nbr=6;
	p = polyfit(log(Dhair(l-nbr:l))-ldh0,log(ppl(l-nbr:l))-leh0,1);
    plot(log(Dhair(l-nbr:l))-ldh0,log(ppl(l-nbr:l))-leh0,'g.');
	powercoeff = [powercoeff p(1)];
    powercoeff2 = [powercoeff2 p(2)];
    %calculons l'erreur
    tmpxx = log(Dhair(l-nbr:l))-ldh0;
    tmpyy = log(ppl(l-nbr:l))-leh0;
    fx=@(x) p(1)*x+p(2);
    eee=sum((tmpyy - fx(tmpxx)).^2 );
    clear tmpyy tmpxx fx
    errcoeff = [errcoeff eee*10];

	plot(log(Dhair)-ldh0,p(1)*log(Dhair)+p(2),'g--');
	%plot(log(Dhair)-ldh0,-2*log(Dhair)-20,'k--');
	%isd = 1 ./ (Dhair .^2)
	%semilogy(Dhair,isd,'.k');
	figure(1);
	clf;
    hold on;
	errorbar([1:length(powercoeff)],powercoeff,errcoeff,'+');
    errorbar([1:length(powercoeff)],powercoeff2,errcoeff,'g.');
	%clear p;
	fprintf('%d/%d',i,length(f));
	%input('next...');
end

%% let's try to do the same for one curve, by changing the sliding value...
figure(1);
clf; hold on;
figure(2)
clf;
% let selec the data we will fit 
irange=[1];

powercoeff=[];
zmean=[];
zstd=[];

n=100;
for i=irange;
    exp=f(i);
    start=1;
    stop = exp.moving_trap.event.appr.stop;
    d=exp.bead_distance(start:stop)-dd(i)+min(dd)+1;
    force= abs(exp.still_trap_force.tangent(start:stop));
    schunk=floor((stop-start)/n);
    X=[];Y=[];Z=[];
	for a=1:15
		m=2*a;
        %m=2*a+13;
		for b=1:30
            
            nbr=d(length(d)-floor((m/2))*schunk)+b*0.2;
            %nbr=b+12;
			%nbr=2*b+18;
			X(a,b)=m;
			Y(a,b)=nbr;
			Ehair=[];
			Dhair=[];
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
				 plot([d(j*schunk) d(j*schunk+lchunk)],[force(j*schunk) force(j*schunk+lchunk)],'ro-');
				 Ehair = [Ehair abs(deltaf*tempdd/deltad)];
				 clear tempdd p;
			end
			
			l = length(Dhair);
			%fitrange=[l-nbr:l];
            fitrange=(Dhair<nbr);
			
			figure(2);
			clf;
			hold on;
			plot(log(Dhair(fitrange)),log(ppl(fitrange)),'bo-');
			plot(log(Dhair),log(ppl),'r-');

			p = polyfit(log(Dhair(fitrange)),log(ppl(fitrange)),1);
			powercoeff = [powercoeff p(1)];
			plot(log(Dhair),p(1)*log(Dhair)+p(2),'g--');
			plot(log(Dhair),-2*log(Dhair)-20,'k--');
			
			figure(1);
			clf;
			plot(powercoeff,'+');
			Z(a,b)=p(1);
			%mesh(X,Y,Z,'EdgeColor','black');
			clear p;
			fprintf('\b%d',b);
			%input('next...');
        end
        %clear powercoeff;
        fprintf('\b\b\b..');
    end
    fprintf('\n---------\n');
    figure(4);
    surf(X,Y,Z);
   % axis([10   25   10   40   -2.5   -1]);
    xlabel('m: average length(raw data)');
    ylabel('nbr: fit length power law');
    zlabel('z: power law exponent');
    colormap hsv;
    colorbar;
    zmean(i)=mean(reshape(Z,1,[]));
    zstd(i)=std(reshape(Z,1,[]));

clear Dhair Ehair a b cps d dd deltad delta dp exp fitrange force irange tm ttm ttm2 


%
%clear X Y Z

s=size(Z);
clear a b tmp stdz Xg Yg
avgz=[];
stdz=[];
Xg=[];
Yg=[];
%[Xg,Yg] = meshgrid(1:1:3);
urange=2;
vrange=2;
for a=[1+urange:s(1)-urange]
	for b=[1+vrange:s(2)-vrange]
        tmp=[];
        for u=-urange:urange
            for v=-vrange:vrange
                tmp = [tmp Z(a+u,b+v)];
            end
        end
        %tmp = [Z(a-1,b-1) Z(a-1,b) Z(a-1,b+1) Z(a,b-1) Z(a,b) Z(a,b+1) Z(a+1,b-1) Z(a+1,b) Z(a+1,b+1)];
		avgz(a-urange,b-vrange)=mean(tmp);
        stdz(a-urange,b-vrange)=std(tmp);
        clear tmp
        Xg(a-urange,b-vrange)=X(a,b);
        Yg(a-urange,b-vrange)=Y(a,b);
	end
end;
clear s
%figure(3);
%clf;
%surf(Xg,Yg,stdz,'EdgeColor','black');
%figure(4);
%clf;
%surf(Xg,Yg,avgz,'EdgeColor','black');
% search least variance extract mean value
avgzr=reshape(avgz,1,[]);
stdzr=reshape(stdz,1,[]);
stdzr_a(i) = min(stdzr);
avgzr_a(i) = avgzr(stdzr==min(stdzr));
ind=find(stdz==min(stdzr));
s=size(stdz);
row = mod(ind,s(1));
col = (ind-row)/s(1)+1;
m_a(i) = X(row+urange,col+vrange);
nbr_a(i) = Y(row+urange,col+vrange);
figure(1)
clf;
plot(nbr_a,m_a,'ro');
xlabel('nbr');
ylabel('m');

%figure(3)
%clf;
%plot(nbr_a,'g*');
%ylabel('nbr');

figure(4)
clf;
errorbar([1:length(avgzr_a)],avgzr_a,stdzr_a,'o');
ylabel('avg');

end



%% tetative rescaling
% on va essayer de traiter toute les courbes de façon à avoir 
% en (0,1) le point avec un maximum de force
% et en (1/2) (1/2) le point avec le maximum de force sur 2

figure(2);
clf
rescaledForce_a = [];
rescaledDistance_a = [];
dd_a = [];
fmax_a =[];
u=length(g);
for j=1:u;
    i=1;
    exp   = g(i);
    start = exp.moving_trap.event.appr.start;
    stop  = exp.moving_trap.event.appr.stop;

    distance = exp.bead_distance(start:stop);
    force    = exp.still_trap_force.tangent(start:stop);

    %determination de la force maximal
    [maxforce,maxindice] = max(force);
    maxdistance          = distance(maxindice);
    fmax_a = [fmax_a maxforce];
    %determination de la force à la moitiée et de la distance correspondante
    frac=1/2;
    demimax = maxforce*frac;

    [~,demiindice] = min(abs(force-demimax));
    demiforce      = force(demiindice);
    demidistance   = distance(demiindice);

    %figure(1);
    %clf
    %hold on;
    %plot(distance,force,'g+');
    %plot([maxdistance demidistance],[maxforce demiforce],'ro');

    figure(2);
    hold on;
    dd_a = [dd_a (demidistance-maxdistance)];
    rescaledDistance = (distance-maxdistance)/(demidistance-maxdistance);
    rescaledForce    = (force/maxforce);
    %plot(rescaledDistance,rescaledForce,':');
    %plot([0 1],[1 frac],'ro');
    fprintf('\b\b\b\b\b\b\b%03d/%03d',i,length(g));
    rescaledDistance_a = [rescaledDistance_a rescaledDistance];
    rescaledForce_a = [rescaledForce_a rescaledForce];
    g(1)=[];
    clear demidistance demiforce demiindice maxforce maxindice maxdistance rescaledDistance 
    clear demimax i start stop force distance ans rescaledForce exp
end

%à trier
tosort = [rescaledDistance_a;rescaledForce_a];
sorted = sortrows(tosort',1)';
clear tosort rescaledDistance_a rescaledForce_a;

sortedDistance = sorted(1,:);
sortedForce = sorted(2,:);

%let's do packet of...
n=100;
dmax=max(sortedDistance(sortedDistance < 50));
step = dmax/n;
for i=1:n
  keep = sortedDistance > (i-1)*step & sortedDistance < i*step ;
  stat.d(i) = mean(sortedDistance(keep));
  stat.f(i) = mean(   sortedForce(keep));
  stat.dstd(i) = std(sortedDistance(keep));
  stat.fstd(i) = std(   sortedForce(keep));
end

%clear sorted sortedDistance step keep
%clear sortedForce dmax n
figure(1)
clf;
hold on;
errorbar(stat.d,stat.f,stat.fstd,'k'); 
plot([0 1],[1 frac],'ro');

%%
figure(3);
range=2:20;
ofset = 1;
%hold on;
p=polyfit(log(stat.d(range)+ofset),log(stat.f(range)),1);
%clf
%plot(log(stat.d(range)),log(stat.f(range)),'ro');


%plot(log(stat.d),log(stat.f),'*');

%
loglog((stat.d+ofset),(stat.f),'o');
hold on
%axis equal
loglog(stat.d+ofset,exp(polyval(p,log(stat.d+ofset))),'k--');












 

