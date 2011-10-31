%% Cut according to previous fits the d0 and the f0, and the rescale force/f_max, and d/d_min
dl0=[];fl0=[];
dl10=[];fl10=[];
dl30=[];fl30=[];
dl50=[];fl50=[];
dlt=[];flt=[];
for i=1:length(fit_res_fin)
    i
    d=(fit_res_fin(i).d-fit_res_fin(i).d0)./min(fit_res_fin(i).d-fit_res_fin(i).d0);
    f=(fit_res_fin(i).force+fit_res_fin(i).f0*1e-12)./max(fit_res_fin(i).force+fit_res_fin(i).f0*1e-12);
    [dl]=linbin(d,10);
    [fl]=linbin(f,10);
    %hold on
    %loglog(dl,fl,'+')
    if fit_res_fin(i).cp==0
       if max(dl)<50
          dl0=[dl0 dl];
          fl0=[fl0 fl];
          plot(dl,fl);
          i;
       else
           continue;
       end
    elseif fit_res_fin(i).cp==10
       dl10=[dl10 dl];
       fl10=[fl10 fl];            
    elseif fit_res_fin(i).cp==30
        if max(fl)<=1
           dl30=[dl30 dl];
           fl30=[fl30 fl];
        else
            continue;
        end
        
    elseif fit_res_fin(i).cp==50
       dl50=[dl50 dl];
       fl50=[fl50 fl]; 
    end
       dlt=[dlt dl];
       flt=[flt fl];
    
end

