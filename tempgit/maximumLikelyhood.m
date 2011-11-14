function ret=maximumLikelyhood(exp)
    alpha_a = 0.3:0.1:2;
    eb=erasableBuffer;
    start=exp.moving_trap.event.appr.start;
    stop=exp.moving_trap.event.appr.stop;
    jump=300;
    xdata=exp.bead_distance(start:jump:stop);
    ydata = exp.still_trap.force.r(start:jump:stop)*(10e11);
    X0=[0,4,1];
    for i=1:length(alpha_a)
        eb.counter(i,length(alpha_a));
        alpha = alpha_a(i);
        f=@(X) cout(@(u)model(u,X(1),X(2),X(3),alpha),xdata,ydata);
        [A,B,C] = fminsearch(f,X0);
        ret.status(i)=C;
        ret.f0(i)=A(1);
        ret.d0(i)=A(2);
        ret.k(i)=A(3);
        ret.error(i)=B;
        ret.alpha(i)=alpha;
        
        clf
        figure(1)
        plot(xdata,ydata,'+');
        hold on 
        pdata=[A(2)+0.05:0.05:max(xdata)];
        plot(pdata,model(pdata,A(1),A(2),A(3),alpha),'--');
        hold off
        drawnow;
        %input('...');
        
    end
    figure(2);
    clf
    plot([ret.alpha],[ret.error],'+');
end

function ret=model(u,fo,d0,k,alpha)
    ret = arrayfun(@(x)fo + k/(x-d0)^alpha,u);
    %ret = arrayfun(@(x)fo + k*exp(-x/d0),u);
end

function ret=cout(f,xdata,ydata)
    ret=sum( (f(xdata)-ydata).^2);
end

