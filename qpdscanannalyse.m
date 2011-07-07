%%
% dans leprogramme d'exemple, lx ly représentent la longueur des scan
% suivant x et y 
% tr représente le time ratio
%%
%load('test_qpd_rawdata.mat');
%data=data_r(:,1:end);%%
%grid_size=20;
%number_of_point=10;
%%
stp = 2/(q.grid_size-1);

data=q.data_r(:,1:end);
xtrap   = data(1,:);
ytrap   = data(2,:);
aodpower= data(3,:);
qpdx    = data(4,:);
sigmax  = data(5,:);
qpdy    = data(6,:);
sigmay  = data(7,:);
qpdsum  = data(8,:);
sigmaqpdsum          = data(9,:);
powerphotodiode      = data(10,:);
sigmapowerphotodiode = data(11,:);

%dx = qpdx ./ qpdsum;
%dy = qpdy ./ qpdsum;


%mx  = mean(reshape(dx,number_of_point,[]));
%my  = mean(reshape(dy,number_of_point,[]));
ms  = mean(reshape(qpdsum,q.number_of_point,[]));
rms = ms/max(ms);
mpd = mean(reshape(powerphotodiode,q.number_of_point,[]));
inpower = mpd /max(mpd);

[X Y] = meshgrid(-1:stp:1,-1:stp:1);

% figure(1)
% surf(X,Y,reshape(mx,grid_size,[]))
% title('mx qpd fonction de la position' )
% 
% figure(2)
% quiver(X,Y,reshape(mx,grid_size,[]),reshape(my,grid_size,[]))
% title('drift photodiode fct de la position');
% 
% figure(3)
% surf(X,Y,reshape(my,grid_size,[]))
% title('my qpd vs position');
% 
% figure(4)
% surf(X,Y,reshape(ms,grid_size,[]))
% title('évolution de la puissance sur la qpd en fonction de la position');
% 
% figure(5)
% surf(X,Y,reshape(inpower,grid_size,[]))
% title('évolution de la puissance sur la PHOTODIODE en %');

figure(11)
Z=reshape(inpower,q.grid_size,[]);
%iz=interp2(Z);
%ix=interp2(X);
%iy=interp2(Y);
[C,h] = contour(X,Y,Z,[0.8:0.02:1.0]);
set(h,'ShowText','on','TextStep',get(h,'LevelStep')*10)
%colormap cool

figure(12)
ZZ=reshape(rms,q.grid_size,[]);
%[C,h] = surf(X,Y,ZZ./Z,[0.8:0.02:1.0]);
surf(X,Y,ZZ./Z)%,[0.8:0.02:1.0]);
%set(h,'ShowText','on','TextStep',get(h,'LevelStep'))

%clear mx my dx dy qpdx qpdy qpdsum

