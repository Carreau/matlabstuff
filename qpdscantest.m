% copier coller et modifié des vi pour voir ce que ça fait 
% notes from connections
% o is from off_time t beginning (10)

% tr is time ratio between trap and scan
% d is distance between scsan point
% scanwidth seem to be an array  (+/-) i x, y, z
% same for cal which is pix to aod
% in is tha trap array (255) mean trap is active

% /!\ don't know why but aom_s is not used
%%
clear
o=10;
tr=10;
d=5;
scan_width =[1 1 0];
cal=[380.0000 257.0000 0.0037 -0.0037 0.0265 -0.0268];
in=[366 283 255 0 0];

%
% on cherche le piège actif
% pos est donc le numero du piège actif
pos = find(in(:,3)==255);

% on a donc la posion initiale du piège actif
x_o = cal(3)*(in(pos,1)-cal(1));
y_o=  cal(4)*(in(pos,2)-cal(2));


%ici on a le voltage du piège actif
aom_t = in(pos,4);
clear in pos

%ici on calcul le déplacement du piège 
dx = cal(5)*d/1000;
dy = cal(6)*d/1000;
clear cal

% on calcul le nombre du point du scan dans chaqu direction
    
n_st_x = round(scan_width(1)*1000/d);
n_st_y = round(scan_width(2)*1000/d);
clear scan_width d

% on calcul les futurs positions du pièges
x_s = x_o+[-n_st_x:n_st_x]*dx;
y_s = y_o+[-n_st_y:n_st_y]*dy;
clear dx dy n_st_x n_st_y

% on les intercales de 'zeros' en concatainant verticalemetn, puis reshape
x_t(1,:)=x_s;
x_t(2:tr+1,:) = x_o*ones(tr,length(x_s));
y_t(1,:)=y_s;
y_t(2:tr+1,:) = y_o*ones(tr,length(y_s));

clear x_s y_s tr

x_i = reshape(x_t,1,size(x_t,1)*size(x_t,2));
y_i = reshape(y_t,1,size(y_t,1)*size(y_t,2));
clear x_t y_t

% we add 'zeros' at the beggining and the end
x_f = [x_i ones(1,length(y_i))*x_o];
y_f = [ones(1,length(x_i))*y_o y_i];
clear x_i y_i

% and ofset it
x_f = [ ones(1,o)*x_o x_f];
y_f = [ ones(1,o)*y_o y_f];

clear x_o y_o o
% create an array with the power.
a_f = aom_t*ones(1,length(x_f));

%



[x_s,y_s] = meshgrid(-1:0.5:1,-1:0.5:1);
x_s=reshape(x_s,1,[]);
y_s=reshape(y_s,1,[]);
x_f=reshape(ones(10,1)*x_s,1,[]);
y_f=reshape(ones(10,1)*y_s,1,[]);
clear x_s y_s;
a_f = aom_t*ones(1,length(x_f));
clear aom_t
