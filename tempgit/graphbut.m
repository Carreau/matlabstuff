%function graphbut(f);
f = load('butdata.csv');
prizelist=unique(sort(f(:,2)));
t0 = min(f(:,1));
f(:,1)=f(:,1)-t0;
figure(1);
clf
hold on 
for itarget=1:length(prizelist)
      target = prizelist(itarget)
      rangs=(f(:,2)==target);
      f(rangs,1)
      f(rangs,3)
      %plot(f(rangs,1),f(rangs,3),'+');
      
end

datelist = unique(sort(f(:,1)));

for itarget=1:length(datelist)
      target = datelist(itarget);
      rangs=(f(:,1)==target);
      mean(f(rangs,3))
      errorbar(itarget+10,mean(f(rangs,3)),std(f(rangs,3)),'ro');
end
tm =[];
mm = [];
for itarget=1:length(datelist)
      target = datelist(itarget);
      rangs=(f(:,1)==target);
      mean(f(rangs,3))
      tm=[tm itarget+10];
      mm=[mm median(f(rangs,3))];
end

plot(tm,mm,'go--');