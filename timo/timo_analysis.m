%This program will use MAtthias' stuff to load the data, and prefilter it.
%Then I will see what I can do with it.
%% load data
g =    temptest('30_mars_10cp_25apr.mat');
%g = [g temptest('30_mars_30cp_25apr.mat')];
%g = [g temptest('30_mars_50cp_25apr.mat')];
%g = [g temptest('8mars25Arp50cp.mat')];
%g = [g temptest('8mars25arp00cp.mat')];
%g = [g temptest('run1_1-mars-2011_25arp-10cp.mat')];
%g = [g temptest('run3_1-mars-2011_25arp-30cp.mat')];
%g = [g temptest('run4_1-mars-2011_25arp-30cp.mat')];



%% demande si garder ou pas la courbe (a la main)
garde=[];
for i=1:length(g)
  % err(i);
   g(i).showfitTau;
   inp=input('garder ? [y/N] ','s');
   if(inp=='n')
	   garde=[garde 0];
	   disp('on garde pas');
   else
	   garde=[garde 1];
	   disp('on garde');
   end
end
%% tetative rescaling
% on va essayer de traiter toute les courbes de façon à avoir 
% en (0,1) le point avec un maximum de force
% et en (1/2) (1/2) le point avec le maximum de force sur 2

figure(2);
clf
rescaledForce_a = [];
rescaledDistance_a = [];
dd_a = [];
fmax_a =[];
for i=1:length(g);
    exp   = g(i);
    start = exp.moving_trap.event.appr.start;
    stop  = exp.moving_trap.event.appr.stop;

    distance = exp.bead_distance(start:stop);
    force    = exp.still_trap_force.tangent(start:stop);

    %determination de la force maximal
    [maxforce,maxindice] = max(force);
    maxdistance          = distance(maxindice);
    fmax_a = [fmax_a maxforce];
    %determination de la force à la moitiée et de la distance correspondante
    frac=1/2;
    demimax = maxforce*frac;

    [~,demiindice] = min(abs(force-demimax));
    demiforce      = force(demiindice);
    demidistance   = distance(demiindice);

    %figure(1);
    %clf
    %hold on;
    %plot(distance,force,'g+');
    %plot([maxdistance demidistance],[maxforce demiforce],'ro');

    figure(2);
    hold on;
    dd_a = [dd_a (demidistance-maxdistance)];
    rescaledDistance = (distance-maxdistance)/(demidistance-maxdistance);
    rescaledForce    = (force/maxforce);
    %plot(rescaledDistance,rescaledForce,':');
    %plot([0 1],[1 frac],'ro');
    fprintf('\b\b\b\b\b\b\b%03d/%03d',i,length(g));
    rescaledDistance_a = [rescaledDistance_a rescaledDistance];
    rescaledForce_a = [rescaledForce_a rescaledForce];
    
    clear demidistance demiforce demiindice maxforce maxindice maxdistance rescaledDistance 
    clear demimax i start stop force distance ans rescaledForce exp
end

%à trier
tosort = [rescaledDistance_a;rescaledForce_a];
sorted = sortrows(tosort',1)';
clear tosort rescaledDistance_a rescaledForce_a;

sortedDistance = sorted(1,:);
sortedForce = sorted(2,:);

%let's do packet of...
n=70;
dmax=max(sortedDistance(sortedDistance < 50));
step = dmax/n;
for i=1:n
  keep = sortedDistance > (i-1)*step & sortedDistance < i*step ;
  stat.d(i) = mean(sortedDistance(keep));
  stat.f(i) = mean(   sortedForce(keep));
  stat.dstd(i) = std(sortedDistance(keep));
  stat.fstd(i) = std(   sortedForce(keep));
end

clear sorted sortedDistance step keep
clear sortedForce dmax n
figure(1)
clf;
hold on;
errorbar(stat.d,stat.f,stat.fstd,'k'); 
plot([0 1],[1 frac],'ro');