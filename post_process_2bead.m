function save_data=post_process_2bead(f_path)

%r is the bead radius, unfortunately, I do not save that one from the
%preogram
r=4.5e-6/2
%this is hte forca over which I try to align all the datasets
f_align=10e-12

%I will give it the preprocesses save_data structure (either the raw, or
%the prefiltere (1µs)
%Then we will seperate the data in the 3 phases (approach, rest,
%retract)(this will be done by just getting the relevant positions).
%Here I will also get the first contact point to be defined by a 5 pN Force
%We will see if we can extract a certain amount of quantitative data:
%1. the dissipated energy (the area under the approach retract
%2. The indentation (get the point where the force was higher than 5 pn, to
%align
%3. The distance dependent elasticity (not taking the viscosity into
%account
%4. The time constant of the exponential decay. And the relative decay with
%reference to indentation
%5. Overlay all the approach curves, with the 'contact point at 0'
%6. Overlay all the retract curves, with the 'contact point at 0'

 if (nargin == 0)
     [file,spath]=uigetfile('E:\Science\data\collection\elasticity');
    f_path=[spath,filesep,file];
 end 
% %load the file
 load(f_path)
 
for i=1:length(save_data)
    %first I calculate the dissipated energy
    save_data(i).diss_energy=1e-6*sum(diff(save_data(i).d)'.*squeeze(save_data(i).f(1,1,2:end)));
    %now get teh points where the approach stated and ended, and where the
    %retraction started
    tr=save_data(i).trap_pos;
    tr_m=sqrt(squeeze(tr(2,1,:)).^2+squeeze(tr(2,2,:)).^2);
    tr_d_0=find(diff(tr_m)==0);
    tr_d_no0=find(diff(tr_m)~=0);
    
    save_data(i).appr_start=tr_d_no0(1);
    save_data(i).retr_start=tr_d_0(end);
    tr_d_0_cut=tr_d_0;
    tr_d_0_cut(find(tr_d_0_cut<=tr_d_no0(1)))=[];
    save_data(i).appr_stop=tr_d_0_cut(1);
    
    %so now I need to estimate the 'touching point'. I will try to do this
    %by getting the point where the forces are bigger than f_align
    [fm,pm]=min(abs(f_align-squeeze(save_data(i).f(1,1,1:save_data(i).appr_stop))));
    save_data(i).touch_point=pm;
   
    %I will also try to get the touching point for the retraction phase 'touching point'. I will try to do this
    %by getting the point where the forces are bigger than f_align
    [fm,pm_r]=min(abs(f_align-squeeze(save_data(i).f(1,1,save_data(i).retr_start:end))));
    pm_r=pm_r+save_data(i).retr_start;
    save_data(i).touch_point_retr=pm_r;
    
    
    %This point is now used to calculate the indentation
    save_data(i).indent=save_data(i).d(pm)-save_data(i).d(save_data(i).retr_start);
    
    %Now we can also use the 'touching point as reference for the elasticy
    %estimate as a function of the indentation. The formular is:
    %F=4/3*E/(1-(nu)^2).*d.^(3/2)*sqrt(r)+f0;
    f_ind=squeeze(save_data(i).f(1,1,pm:save_data(i).appr_stop))';
    d_ind=save_data(i).d(pm:save_data(i).appr_stop);
    d_ind=-(d_ind-d_ind(1))*1e-6;
    save_data(i).E=(f_ind-f_align)*3*(1-0.5^2)./(4*d_ind.^(3/2)*sqrt(r));
    save_data(i).E_final=mean(save_data(i).E(round(numel(save_data(i).E)/2:end)));
    
    %Next is the exponential decay first extract the regime from the data,
    %then fit an exponential of the form d(t)=alpha*exp(-t/tau)+d_e
    
    d_t=save_data(i).d(save_data(i).appr_stop:save_data(i).retr_start);
    t_t=[0:length(d_t)-1]*1/save_data(i).parameters.Effective_Sampling_Rate.value;
    
    start_point = [d_t(1)-d_t(end),d_t(end),t_t(end)];
    %options=optimset('Display','iter');
    %estimates = fminsearch(@exp_dec, start_point,options,t_t,d_t);
    estimates = fminsearch(@exp_dec, start_point,[],t_t,d_t);
    [sse,exp_fit]=exp_dec(estimates,t_t,d_t);
    estimates;
    save_data(i).alpha=estimates(1);
    save_data(i).tau=estimates(3);
    save_data(i).d0=estimates(2);
    save_data(i).rel_decay=estimates(1)/estimates(2);
    
    %finally I will get the right regimes for the overlay replotting
    save_data(i).f_app_overlay=squeeze(save_data(i).f(1,1,1:save_data(i).appr_stop));
    save_data(i).f_retr_overlay=squeeze(save_data(i).f(1,1,save_data(i).retr_start:end));
    save_data(i).d_app_overlay=save_data(i).d(1:save_data(i).appr_stop)-save_data(i).d(pm);
    save_data(i).d_retr_overlay=save_data(i).d(save_data(i).retr_start:end)-save_data(i).d(pm_r);

end
sa_path=[f_path(1:end-4),'_postprocessed.mat'];

%and finally wa save in the same file with 'post_prossed' prepend
save(sa_path,'save_data');

function [sse, FittedCurve] =exp_dec(params,t,d)
        %global eta2 fc
        a = params(1);
        b = params(2);
        tau= params(3);
        FittedCurve=a*exp(-t./tau)+b;
        ErrorVector = (FittedCurve - d)./d;
        ErrorVector(isnan(ErrorVector))=[];
        sse = abs(sum(abs(ErrorVector) .^ 2))*1e10;        
        
