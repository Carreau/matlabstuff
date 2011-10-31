function get_power_fit_on_experiment_2()
 
str='alldata-y.mat';
eb = erasableBuffer;
f =load(str);
m=f.v
	for i=1:3;
        eb.counter(i,length(m));
		exp=m(i);
		start= exp.start;
		stop = exp.stopt;
		d=exp.d(start:stop);
		force= abs(exp.f(start:stop));
		
		% in this version, I estimate the rage that should be used for the fit,
		% if the data is steep, only use 10 times the x range after the data
		% has fallen to 20% of the max value
		fmax=max(force);
		f_end=mean(force(1:ceil(end*.01)));
		f_lim=(fmax-f_end)*.2;
		[~,p]=min(abs(force-f_end-f_lim));
		%find correspondinf d, and check if the 20 time the d is smaller than
		%the final d
		dp=d(p);
		d_min=min(d);
		de=10*(dp-d_min)+d_min;
		if de>max(d)
			p_fin=1;
		else
			[~,p_fin]=min(abs(d-de));
		end
	
		
		[d0(i),f0(i),alpha(i),k(i),f_out,f_gof]=fit_power_with_offsets(d(p_fin:end),force(p_fin:end)*1e12);%#ok<AGROW>
		fit_res(i).fit=f_out; %#ok<AGROW>
		fit_res(i).gof=f_gof;%#ok<AGROW>
		fit_res(i).d0=d0(i);%#ok<AGROW>
		fit_res(i).f0=f0(i);%#ok<AGROW>
		fit_res(i).alpha=alpha(i);%#ok<AGROW>
		fit_res(i).k=k(i);%#ok<AGROW>
		fit_res(i).d=d;%#ok<AGROW>
		fit_res(i).force=force;%#ok<AGROW>
		%fit_res(i).cp=exp.cp;%#ok<AGROW>
		%fit_res(i).arp=exp.arp;%#ok<AGROW>
        %fit_res(i).time_m=exp.time_m;%#ok<AGROW>
	end
		
	 %and now save the stuff
	  
save([str(1:end-4),'_fit_res.mat'],'fit_res')
system('say matlab script terminated without error');
end