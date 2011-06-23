%%
%f= load('calibration_simons1000_AOM=0_bead=1.mat');
figure(1);clf
figure(2);clf
dat.x= [];
dat.y=[];
cd ../calibration1001/

f= load('calibration_simons1000_AOM=0_bead=1.mat');
dat.x= f.x-mean(f.x);
dat.y= f.y-mean(f.y);

f= load('calibration_simons1001_AOM=0_bead=1.mat');
dat.x = [dat.x [f.x-mean(f.x)] ];
dat.y = [dat.y [f.y-mean(f.y)] ];

f= load('calibration_simons1002_AOM=0_bead=1.mat');
dat.x = [dat.x [f.x-mean(f.x)] ];
dat.y = [dat.y [f.y-mean(f.y)] ];

cd ../calibration1002/

f= load('calibration_simons1000_AOM=0_bead=1.mat');
dat.x = [dat.x [f.x-mean(f.x)] ];
dat.y = [dat.y [f.y-mean(f.y)] ];

f= load('calibration_simons1001_AOM=0_bead=1.mat');
dat.x = [dat.x [f.x-mean(f.x)] ];
dat.y = [dat.y [f.y-mean(f.y)] ];

f= load('calibration_simons1002_AOM=0_bead=1.mat');
dat.x = [dat.x [f.x-mean(f.x)] ];
dat.y = [dat.y [f.y-mean(f.y)] ];

cd ../calibration1003/

f= load('calibration_simons1000_AOM=0_bead=1.mat');
dat.x = [dat.x [f.x-mean(f.x)] ];
dat.y = [dat.y [f.y-mean(f.y)] ];

f= load('calibration_simons1001_AOM=0_bead=1.mat');
dat.x = [dat.x [f.x-mean(f.x)] ];
dat.y = [dat.y [f.y-mean(f.y)] ];

f= load('calibration_simons1002_AOM=0_bead=1.mat');
dat.x = [dat.x [f.x-mean(f.x)] ];
dat.y = [dat.y [f.y-mean(f.y)] ];

%
d = reshape([[dat.x] [dat.y]],length(dat.y),2);

%
figure(1)
hist3(d,[70 70]);
n = hist3(d,[100 100]);
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');

% admettonsn
figure(2)
n1 = n'; 
n1( size(n,1) + 1 ,size(n,2) + 1 ) = 0;
contour(n1,15)
set(gca, 'DataAspectRatio', [1 1 1])