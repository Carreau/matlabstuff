function get_power_fit_on_experiment()
%% 
%dat(1).name='8mars25Arp50cp.mat';
dat(1).name='8mars25arp00cp.mat';
%dat(3).name='run1_1-mars-2011_25arp-10cp.mat';
%dat(4).name='run3_1-mars-2011_25arp-30cp.mat';
%dat(5).name='run4_1-mars-2011_25arp-30cp.mat';
%dat(6).name='30_mars_10cp_25apr.mat';
%dat(7).name='30_mars_30cp_25apr.mat';
%dat(8).name='30_mars_50cp_25apr.mat';

%%
eb = erasableBuffer;
for ji=1:length(dat)
    disp('loop');
	clear f g d force exp f_out
	%analyze the bead approach data
	str=[dat(ji).name];
	g =    temptest(str);
	%g = [g temptest('30_mars_30cp_25apr.mat')];
	%g = [g temptest('30_mars_50cp_25apr.mat')];
	%g = [g temptest('8mars25Arp50cp.mat')];
	%g = [g temptest('8mars25arp00cp.mat')];
	%g = [g temptest('run1_1-mars-2011_25arp-10cp.mat')];
	%g = [g temptest('run3_1-mars-2011_25arp-30cp.mat')];
	%g = [g temptest('run4_1-mars-2011_25arp-30cp.mat')];
	% check if sum of force is constent during the resting
	clear i s m somme
	%s.moy.x=[];
	%m.moy.x=[];
	%s.ecr.x=[];
	%m.ecr.x=[];
	%s.moy.y=[];
	%m.moy.y=[];
	%s.ecr.y=[];
	%m.ecr.y=[];

	%somme.moy.x=[];
	%somme.std.x=[];
	%somme.moy.y=[];
	%somme.std.y=[];
	err=[];
	nn=[];
	clear f;
	tic;
    %somme=zeros(length(g));
    %s=zeros(length(g));
    %m=zeros(length(g));
    %fit_res=zeros(length(g));
	for i=1:length(g)
        eb.counter(i,length(g));
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
		
		[~,~]=hist(sty-mean(sty),n);
		%plot(b,a,'r.');
		%hold on
		[~,~]=hist(mty-mean(mty),n);
		%plot(b,a,'b.');
		[a,b]=hist(sty+mty-mean(sty+mty),n);
		%[mu,sig,ermu,ersig] = normfit(sty+mty-mean(sty+mty));

		%plot(b,a,'g*--');
		%fu= @(x)normpdf(x,mu,sig)*max(a)/normpdf(0,0,sig);
		%plot(b,fu(b),'k--','LineWidth',1.2);
		%err= [err sum((a-fu(b)).^2)/length(mty)];
		
		hold off;
		
		s(i).moy.y=mean(sty);
		m(i).moy.y=mean(mty);
		s(i).ecr.y=std(sty);
		m(i).ecr.y=std(mty);
		
		somme(i).moy.x=mean(stx+mtx);
		somme(i).std.x=std(stx+mtx);
		somme(i).moy.y=mean(sty+mty);
		somme(i).std.y=std(sty+mty);
		

		%fprintf('\n%d over %d ',i,length(g));

		hold off;

		%input('next...');
	end
	clear n nn range 
	toc
	x=1:length(g);


	for i=1:length(g);
        eb.counter(i,length(g));
		exp=g(i);
		start=1;
		stop = exp.moving_trap.event.appr.stop;
		d=exp.bead_distance(start:stop);
		force= abs(exp.still_trap_force.tangent(start:stop));
		
		% in this version, I estimate teh rage that should be used for the fit,
		% if the data is steep, only use 10 times the x range after the data
		% has fallen to 20% of the max value
		fmax=max(force);
		f_end=mean(force(1:ceil(end*.01)));
		f_lim=(fmax-f_end)*.2;
		[v,p]=min(abs(force-f_end-f_lim));
		%find correspondinf d, and check if the 20 time the d is smaller than
		%the final d
		dp=d(p);
		d_min=min(d);
		de=10*(dp-d_min)+d_min;
		if de>max(d)
			p_fin=1;
		else
			[v,p_fin]=min(abs(d-de));
		end
		
		
		
		[d0(i),f0(i),alpha(i),k(i),f_out,f_gof]=fit_power_with_offsets(d(p_fin:end),force(p_fin:end)*1e12);
		fit_res(i).fit=f_out;
		fit_res(i).gof=f_gof;
		fit_res(i).d0=d0(i);
		fit_res(i).f0=f0(i);
		fit_res(i).alpha=alpha(i);
		fit_res(i).k=k(i);
		fit_res(i).d=d;
		fit_res(i).force=force;
		fit_res(i).cp=exp.cp;
		fit_res(i).arp=exp.arp;
        fit_res(i).time_m=exp.time_m;
	end
		
	 %and now save the stuff
	  
	save([str(1:end-4),'_fit_res_test.mat'],'fit_res')

end
%system('say matlab script terminated without error');
end