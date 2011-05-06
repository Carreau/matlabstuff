%%
%f= load('calibration_simons1000_AOM=0_bead=1.mat');

%
dat = reshape([[f.x-mean(f.x)] [f.y-mean(f.y)]],length(f.y),2);

%
%hist3(dat,[50 50]);
n = hist3(dat,[50 50]);
%set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');

% admettonsn
n1 = n'; 
n1( size(n,1) + 1 ,size(n,2) + 1 ) = 0;
contour(n1,3)