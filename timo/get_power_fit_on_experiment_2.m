function [fit_res]=get_power_fit_on_experiment_2(struct,relaxed)
% if relaxed set to true, fit wont try to stay at d0 >0
%

%fit_res = 
arrayfun(@(x)corefunc(x,relaxed),struct);
		
	 %and now save the stuff
	 
%save([datestr(now,'YYYY-mm-DD_HH-MM-ss'),sprintf('-relaxed_%d_',relaxed),'_fit_res.mat'],'fit_res')
system('say matlab script terminated without error');
end

function corefunc(ex,relaxed)
        %let's start parralisation crazyness
        %get the uuid and construc a filename from it
        uuid     = char(ex.UUID.toString());
        basename=sprintf('parralelisme_crazyness_%s_relaxed_%d',uuid,relaxed);
        lockfile = sprintf('%s.lock',basename); 
        matfile = sprintf('%s.mat',basename); 
        if (exist(matfile,'file'))
           disp('matfile exist')
           return 
        end
        if (exist(lockfile,'file'))
           disp('lockfile exist, skipping')
           return 
        end
        
        save(lockfile,'-ascii','uuid')
        
        persistent niter;
        %niter=niter+1;
        if isempty(niter)
           disp('is empty');
           niter=0;
        else
             niter=niter+1;
        end
        fprintf('iteration %d\n', niter)
		start= ex.start;
		stop = ex.stopt;
		d=ex.d(start:stop);
		
        %take the absolute value, otherwithe some part of the fi crash
        force = abs(ex.f(start:stop));
		
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
		[d0,f0,alpha,k,f_out,f_gof]=fit_power_with_offsets(d(p_fin:end),force(p_fin:end)*1e12,relaxed);
		fit_res.fit = f_out; 
		fit_res.gof = f_gof;
		fit_res.d0  = d0;
		fit_res.f0  = f0;
		fit_res.alpha = alpha;
		fit_res.k   = k;
		fit_res.d   = d;
		fit_res.force = force;
        fit_res.raw = ex;
        fit_res.relaxed = relaxed;
        fit_res.uuid    = char(ex.UUID.toString());
        save(matfile,'fit_res')
        delete(lockfile)
		%fit_res.cp=exp.cp;%#ok<AGROW>
		%fit_res.arp=exp.arp;%#ok<AGROW>
        %fit_res.time_m=exp.time_m;%#ok<AGROW>
end