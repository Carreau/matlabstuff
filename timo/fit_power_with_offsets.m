function [d0,alpha,k,f_out,f_gof]=fit_power_with_offsets(d,f,relaxed)

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

start_point = min(d);%,min(f)];
end_point = max(d);

%can't we just use modelfun ? 
model = @(p)modelfun(p);
if(relaxed)
    lowbound = -10*end_point;
else 
    lowbound = 0;
end
%hi_bound= (start_point+end_point*9)/10
hi_bound= start_point;
estimates = fminbnd(model, lowbound,hi_bound );

d0 = estimates(1);
%f0 = estimates(2)
%eta = estimates(2)


[sse,alpha,k,f_out,f_gof] = model(estimates);
plot(xdata,ydata);
hold on
plot(xdata, k*(xdata-d0).^alpha,'r');
hold off

%why do we pause ?




%-----------------------------------------------------------------------
 function [sse,  alpha,k,f_out,f_gof] = modelfun(params)
     %global sse_old
        persistent sse_old

        %d0 = abs(params(1));
        d0 = params(1);
        %d0 = min(d0,max(xdata')*0.9+min(xdata')*0.1);
        
        posx=find(xdata-d0<=0);
        
        xdata_i=xdata';
        xdata_i(posx)=[];
        ydata_i=ydata';
        ydata_i(posx)=[];
        
        
        %figure(1)
        %fprintf('plotting...for d0:%d versus %d\n',d0,max(xdata')*0.9+min(xdata')*0.1)
        %plot(xdata_i-d0,ydata_i)
        %pause(.1)
        % need to check we have enough datapoint
        if (d0<=min(xdata))
            [f_out,f_gof]=fit((xdata_i-d0), ydata_i,'power1','Lower',[0,-5.001],'Upper',[1e72,-0.0]);
            sse = f_gof.sse;
            alpha=f_out.b;
            k=f_out.a;

        else
            disp('is not  in relaxed mode, cause using sse_old')
            sse = sse_old;
        end

        sse_old=sse;
    end
end