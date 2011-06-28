% copier coller et modifié des vi pour voir ce que ça fait 
% notes from connections
% o is from off_time t beginning (10)

% tr is time ratio between trap and scan
% d is distance between scsan point
% scanwidth seem to be an array  (+/-) i x, y, z
% same for cal which is pix to aod
% in is tha trap array (255) mean trap is active

%%
%clear
%o=10;
%tr=10;
%d=5;
%scan_width =[1 1 0];
%cal=[380.0000 257.0000 0.0037 -0.0037 0.0265 -0.0268];
%in=[366 283 255 0 0];

%%
% on cherche le piège actif
% pos est donc le numero du piège actif
save('D:\temp\matthias\test_qpd_parameters')
pos = find(in(:,3)==255);

%ici on a le voltage du piège actif
aom_t = in(pos,4);
clear in pos


%step size on the grid
grid_size=20;
number_of_point=100;
stp = 2/(grid_size-1);
length([-1:stp:1])

[x_s,y_s] = meshgrid(-1:stp:1,-1:stp:1);
x_s=reshape(x_s,1,[]);
y_s=reshape(y_s,1,[]);
x_f=reshape(ones(number_of_point,1)*x_s,1,[]);
y_f=reshape(ones(number_of_point,1)*y_s,1,[]);
clear x_s y_s;
a_f = aom_t*ones(1,length(x_f));
clear aom_t
