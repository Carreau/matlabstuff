function [xo,yo,lb,lb_std,lb_ste]=average(x,y,dx)

    l(:,1)=x;
    l(:,2)=y;
    ls=sortrows(l,1);
    pa=1;
    %Now I cut into chunks and average
    for i=1:upper(max(x)/dx)
        [v,p]=min(abs(ls(:,1)-min(x)-i*dx));
        li=ls(1:p,:);
        ls(1:p,:)=[];
        %pa=p;
        lb(:,i)=mean(li,1);
        lb_std(:,i)=std(li,1);
        lb_ste(:,i)=std(li,1)./length(li);
    end
    xo=lb(1,:);
    yo=lb(2,:);
end