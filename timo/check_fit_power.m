%First I will generate 10 curves
clear y x
x=[6.5:0.01:30];
al=-1


for i=1:10
    x(i,:)=x(1,:);
    f0(i)=rand*2;
    k(i)=1+rand*10;
    d0(i)=4.5+rand*1.7;
    y(i,:)=f0(i)+k(i)*(x(i,:)-d0(i)).^al+randn(1,length(x))*.3;
end

%% Next we try to recover the stuff

for i=1:10
    [d0_out(i),f0_out(i),a(i)]=fit_power_with_offsets(x(i,:),y(i,:))

    
end

%% Now with all these values, I will try to collapse the courves, so I will
%%use the f0, and d0 to cut the curve, and the normalize as we did before
clear y_n x_n
for i=1:10
    y_n(i,:)=y(i,:)-f0_out(i);
    x_n(i,:)=x(i,:)-d0_out(i);
    
    y_n(i,:)=y_n(i,:)/mean(y_n(i,1:50));
    [v,p]=min(abs(y_n(i,:)-.5));
    
    x_n(i,:)=x_n(i,:)/x_n(i,p);
   % plot(x_n(i,:),y_n(i,:));
    hold on
end
hold off