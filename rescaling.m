%% chargement des donnée utilies
clear
g =    temptest('30_mars_10cp_25apr.mat');
%g = [g temptest('30_mars_30cp_25apr.mat')];
%g = [g temptest('30_mars_50cp_25apr.mat')];
%g = [g temptest('8mars25Arp50cp.mat')];
%g = [g temptest('8mars25arp00cp.mat')];
%g = [g temptest('run1_1-mars-2011_25arp-10cp.mat')];
%g = [g temptest('run3_1-mars-2011_25arp-30cp.mat')];
%g = [g temptest('run4_1-mars-2011_25arp-30cp.mat')]; %#ok<NASGU>

disp('done loading...');
%% polynome de rescaling




%% generation de données arbitraire pour le rescaling.
%  pour se faire on choisit une loi avec n paramètre
%  arbitraires, puis on stocke les valeurs voulues dans  un 'g' correcte. 
%  il faut stocker 
%    start = exp.moving_trap.event.appr.start;
%    stop  = exp.moving_trap.event.appr.stop;
%    distance = exp.bead_distance(start:stop)-ofset;
%    force    = exp.still_trap_force.tangent(start:stop);
% 
% first the parameters of the laws
% number of false sample

% n           = nombre de courbes à génnérer
% klist       = Force proportionnality factor
% dlist       = 1+0.1*rand(n,1);
% dofset      = 1+0.1*rand(n,1);
% stoplist    = 3.0+1.5*rand(n,1);
% startlist   = 55+3*rand(n,1);
% forceoffset = rand(n,1)/500;

n           = 50;
noise       =      (  1/900          );
klist       = rand(n,1);
dlist       = 1  + ( 0.1*rand(n,1)    );
dofset      =      ( 1+0.1*rand(n,1)  );
stoplist    = 3  + ( 1.5*rand(n,1)    );
startlist   = 25 + ( 25+3*rand(n,1)   );
forceoffset =      0.5*( rand(n,1)/500    );
% m==siMulation , t== Theoriquz, f== fit
mexposant   = 2.0;
texposant   = mexposant;


% data point generated every 'step' in distance
step=-0.005;

figure(1)
clf
clc
title('generated curves')
xlabel('distance');
ylabel('force');

hold on;
clear g;
g=[];
for i=1:n
    bd=(startlist(i):step:stoplist(i));
    g(i).moving_trap.event.appr.start=1; %#ok<SAGROW>
    g(i).moving_trap.event.appr.stop=length(bd);%#ok<SAGROW>
    %fprintf('len: %d \n',length(bd)/5000);
    law=@(x) klist(i)/(dlist(i)*x-dofset(i))^mexposant+forceoffset(i);
    g(i).bead_distance=bd+4.5;%#ok<SAGROW>
    g(i).cp=333;%#ok<SAGROW>
    g(i).still_trap_force.tangent = arrayfun(law,bd)+randn(1,length(bd))*noise;%#ok<SAGROW>
    figure(1)
    plot(g(i).bead_distance(1:10:end),g(i).still_trap_force.tangent(1:10:end));
end



%% tetative rescaling
% on va essayer de traiter toute les courbes de façon à avoir 
% en (0,1) le point avec un maximum de force
% et en (1/2) (1/2) le point avec le maximum de force sur 2
fexposant   = 1.0;
texposant   = fexposant;

%figure(2);
u=length(g);
clf
rescaledForce_a = [];
rescaledDistance_a = [];
dd_a    = zeros(u,1);
dmax_a  = zeros(u,1);
dfrac_a = zeros(u,1);
dquart_a= zeros(u,1);
fmax_a  = zeros(u,1);
cp_a    = zeros(u,1);

%


%ofset du à l'épaiseur de la bille
ofset=4.5;
sstp=1;
fprintf('\n000/%03d:(000/000)',u);
keeped=0;
throwed=0;
interv=0.9:0.1:7;
figure(2);
clf
hold on;
frac=1/3;
upperlim = @(u) arrayfun(@(x) min(0.4./ (0.8*x-0.5)+0.3,1.2),u);
lowerlim = @(u) arrayfun(@(x)   1./ x    -0.2,u);
plot(interv,upperlim(interv),'r--');
plot(interv,lowerlim(interv),'r--');
plot([1 1/frac^(1/fexposant)],[1 frac],'go');
mm=zeros(u,1);
for j=1:u;
    i=j;
    
    exp   = g(i);
    start = exp.moving_trap.event.appr.start;
    stop  = exp.moving_trap.event.appr.stop;

    distance = exp.bead_distance(start:stop)-ofset;
    force    = exp.still_trap_force.tangent(start:stop);
    % on va retirer arbitrairement 4/5 de la moyenne de la force au plus
    % loin
    % force = force- mean(force(start:end/20));
    ee=floor(length(force)/20);
    mm(i)=mean(force(1:ee));
    %fprintf('mean: %d(%d) \n',mm(i),ee-1);
    
    % on adoucis les courbes et on diminue le nombre de points
    distance=linbin(distance,50);
    force=linbin(force,50);
    % distance = slidingavg(distance,100);
    % force = slidingavg(force,100);
        

    %determination de la force maximal
    [fmax,maxindice] = max(force);
    dmax          = distance(maxindice);
    dmax_a(i) = dmax;
    fmax_a(i) = fmax;

    %determination de la force à la moitiée et de la distance correspondante
    
    fracmax = fmax*frac;
    quartmax = fmax*1/4;

    [~,demiindice] = min(abs(force-fracmax));
    [~,quartindice] = min(abs(force-quartmax));
    ffrac      = force(demiindice);
    dfrac   = distance(demiindice);
    quartdistance  = distance(quartindice);
    
    %figure(1);
    %title('target point of current curves')
    %plot(distance(1:sstp:end),force(1:sstp:end),'g+');
    %plot([dmax dfrac],[fmax ffrac],'ro');

    figure(2);
    title('rescaling and filtering');
    hold on;
    %dd_a = [dd_a (demidistance-maxdistance)];
    dfrac_a(i)  = dfrac;
    dquart_a(i) = quartdistance;
    cp_a(i)     = exp.cp;
    
    % on fait le 'rescaling'
    %rescaledDistance = (distance-maxdistance)/(demidistance-maxdistance)+1;
    %a = (1 - 1/frac)/(dmax-dfrac);
    
    a = (1 - frac^(-1/fexposant))/(dmax-dfrac);
    b = (dmax*frac^(-1/fexposant)-dfrac)/(dmax-dfrac);
    resacalingDistancePolynome = [a b];
    clear a b;
    %rescaledDistance = (distance-maxdistance)/(demidistance-maxdistance)+1;
    rescaledDistance = polyval(resacalingDistancePolynome,distance);
    rescaledForce    = (force/fmax);
    
    % on compense la remise à zero precédente
    %rescaledForce = rescaledForce+1/mean(rescaledDistance(1:end/20));
    
    w=rescaledDistance > 20;
    rescaledDistance(w)=[];
    rescaledForce(w)=[];
    % let's remove the curve that really don't match.
    % ie rescales force > 0.3 , au dela de 0.6, où < 0
    w1 = false;%rescaledDistance > 6;
    w2 = false;%rescaledForce    >0.3;
    w3 = false;%rescaledForce    <-0.1;
    w4 = rescaledForce > upperlim(rescaledDistance);
    w5 = rescaledForce < lowerlim(rescaledDistance);
    throw = false; (w1 & (w2 | w3) | w4 |w5);
    if (sum(throw)<1)
        keeped=keeped+1;
        plot(rescaledDistance(1:sstp:end),rescaledForce(1:sstp:end),'-');
        rescaledDistance_a = [rescaledDistance_a rescaledDistance]; %#ok<AGROW>
        rescaledForce_a = [rescaledForce_a rescaledForce];%#ok<AGROW>
    else
        throwed=throwed+1;
        plot(rescaledDistance(1:sstp:end),rescaledForce(1:sstp:end),'r-');
    end
    
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%03d/%03d:(%03d/%03d)',i,length(g),keeped,throwed);
    %g(1)=[];
    clear demidistance demiforce demiindice maxforce maxindice maxdistance rescaledDistance 
    clear demimax i start stop force distance ans rescaledForce exp
    
end
clear interv j keep keeped n ofset p quartdistance quartindice 
clear thro* u w* sstp quartmax

% à trier
tosort = [rescaledDistance_a;rescaledForce_a];
sorted = sortrows(tosort',1)';
clear tosort rescaledDistance_a rescaledForce_a;

sortedDistance = sorted(1,:);
sortedForce = sorted(2,:);

%let's do packet of...
n=50;
clear stat
dmax=max(sortedDistance(sortedDistance < 50));
step = dmax/n;
for i=1:n
  keep = sortedDistance > (i-1)*step & sortedDistance < i*step ;
  stat.d(i) = mean(sortedDistance(keep));
  stat.f(i) = mean(   sortedForce(keep));
  stat.dstd(i) = std(sortedDistance(keep));
  stat.fstd(i) = std(   sortedForce(keep));
  clear keep
end

%clear sorted sortedDistance step keep
%clear sortedForce dmax n
figure(2)
%clf;
hold on;
errorbar(stat.d,stat.f,stat.fstd,'k'); 
plot([1 1/frac^(1/fexposant)],[1 frac],'ro');

%
figure(3);
clf
range=1:30;
ofset = 0;
%hold on;
x=log(stat.d(range)+ofset);
y=log(stat.f(range));

w=not(isnan(x));
range=w;
p=polyfit(x(w),y(w),1);
fprintf('coeff dirrecteur : %d',p(1));
clf
%%%loglog(stat.d(range),stat.f(range),'r+');


%plot(log(stat.d),log(stat.f),'*');

%
loglog((stat.d+ofset),(stat.f),'o');
hold on 
loglog((stat.d+ofset),(stat.f-stat.fstd),'--');
loglog((stat.d+ofset),(stat.f+stat.fstd),'--');
hold on
%axis equal
clear exp
%%%loglog(stat.d+ofset,exp(polyval(p,log(stat.d+ofset))),'k--');
%
r=frac;
xx=1:0.1:15;
for exposant=0.5:0.5:3
    ff=@(x) arrayfun(@(u)1/(u^exposant),x);
    %h=@(x) arrayfun(@(u)1/(u^(1/exposant)),x);
    rescales=@(u) arrayfun(@(x) (x-1)*(1-r^(-1/exposant))/(1-r^(-1/fexposant))+1,u);
    loglog(xx,ff(rescales(xx)),'g--');
end

hold on 

%

f=@(x) arrayfun(@(u)1/(u^texposant),x);
h=@(x) arrayfun(@(u)1/(u^(1/texposant)),x);
%rescales=@(u) arrayfun(@(x) (x-1)*(r^(-1/fexposant)-1)/(r^(-1/texposant)-1)+1,u);
%rescales=@(u) arrayfun(@(x) (x*(1-r^(-1/fexposant))+(r^(-1/fexposant)-r^(-1/texposant))) / (1-r^(-1/texposant)),u);
rescales=@(u) arrayfun(@(x) (x-1)*(1-r^(-1/texposant))/(1-r^(-1/fexposant))+1,u);

loglog(xx,f(rescales(xx)),'r--');

loglog([1 1/frac^(1/fexposant)],[1 frac],'kx');
xlabel('log distance (rescaled)');
ylabel('log force (rescaled)');

