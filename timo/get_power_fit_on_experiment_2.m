function [fit_res]=get_power_fit_on_experiment_2(struct,relaxed)
% if relaxed set to true, fit wont try to stay at d0 >0
%
eb = erasableBuffer;
m=struct;
	for i=1:length(m);
        eb.counter(i,length(m));
		exp=m(i);
		start= exp.start;
		stop = exp.stopt;
		d=exp.d(start:stop);
		
        %take the absolute value, otherwithe some part of the fi crash
        force = abs(exp.f(start:stop));
		
		% in this version, I estimate the range that should be used for the fit,
		% if the data is steep, only use 10 times the x range after the data
		% has fallen to 20% of the max value
		fmax = max(force);
        
		f_end = mean(force(1:ceil(end*.01)));
        
		f_lim =(fmax-f_end)*.2;
		
        % p is know the index of the force that is the closer from the 20%
        % max
        [~,p]=min(abs(force-f_end-f_lim));
        
		% find correspondinf d, and check if the 20 time the d is smaller
		% than the final d
        % (Matt: not sure I understand...) 
		dp=d(p);
		d_min=min(d);
		de=10*(dp-d_min)+d_min;
		if de>max(d)
			p_fin=1;
		else
			[~,p_fin]=min(abs(d-de));
		end
	
		%here we fit the power law
		[d0(i),f0(i),alpha(i),k(i),f_out,f_gof]=fit_power_with_offsets(d(p_fin:end),force(p_fin:end)*1e12,relaxed);%#ok<AGROW>
		fit_res(i).fit=f_out; %#ok<AGROW>
		fit_res(i).gof=f_gof;%#ok<AGROW>
		fit_res(i).d0=d0(i);%#ok<AGROW>
		fit_res(i).f0=f0(i);%#ok<AGROW>
		fit_res(i).alpha=alpha(i);%#ok<AGROW>
		fit_res(i).k=k(i);%#ok<AGROW>
		fit_res(i).d=d;%#ok<AGROW>
		fit_res(i).force=force;%#ok<AGROW>
        fit_res(i).raw=m(i);%#ok<AGROW>
        fit_res(i).relaxed = relaxed;%#ok<AGROW>
        fit_res(i).uuid    = m(i).UUID.toString();%#ok<AGROW>
		%fit_res(i).cp=exp.cp;%#ok<AGROW>
		%fit_res(i).arp=exp.arp;%#ok<AGROW>
        %fit_res(i).time_m=exp.time_m;%#ok<AGROW>
	end
		
	 %and now save the stuff
	 
save([datestr(now,'YYYY-mm-DD_HH-MM-ss'),sprintf('-relaxed_%d_',relaxed),'_fit_res.mat'],'fit_res')
system('say matlab script terminated without error');
end