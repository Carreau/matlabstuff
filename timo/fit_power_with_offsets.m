function [d0,f0,alpha,k,f_out,f_gof]=fit_power_with_offsets(d,f,relaxed)

% might want to use a persistant variable instead of a global. 
% persistant will have the advantage of not varing between calls, but will
% only be visible from inside the scope of the fonction. 
%global sse_old

%I will look at the loglog +d0,f0 and then try to do a linear fit. I will
%then vary the d0,f0 to get the closest to a line as possible.
%if (relaxed == true)
%    disp('fitting with relaxed condition')
%else
%    disp('fitting with unrelaxed condition')
%end

xdata=d;
ydata=f;

start_point = [min(d)];%,min(f)];

%can't we just use modelfun ? 
model = @(p)modelfun(p,relaxed);

estimates = fminsearch(model, start_point);


d0 = estimates(1);
%f0 = estimates(2)
%eta = estimates(2)


[sse,alpha,f0,k,f_out,f_gof] = model(estimates);
plot(xdata,ydata);
hold on
plot(xdata, k*(xdata-d0).^alpha+f0,'r');
hold off

%why do we pause ?
%pause(.5)



%-----------------------------------------------------------------------
 function [sse,  alpha,f0,k,f_out,f_gof] = modelfun(params,relaxed)
     %global sse_old
        persistent sse_old

        d0 = abs(params(1));
        
        posx=find(xdata-d0<=0);
        
        xdata_i=xdata';
        xdata_i(posx)=[];
        ydata_i=ydata';
        ydata_i(posx)=[];
        
        
        if ((d0<=min(xdata) && d0>0 ) || relaxed )
            [f_out,f_gof]=fit((xdata_i-d0), ydata_i,'power2','Lower',[0,-20,-inf],'Upper',[Inf,0,inf]);
            sse = f_gof.sse;
            alpha=f_out.b;
            f0=f_out.c;
            k=f_out.a;

        else
            disp('is not  in relaxed mode, cause using sse_old')
            sse = sse_old;
        end

        sse_old=sse;
    end
end