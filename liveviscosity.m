function [files] =  liveViscosity()
	recursiveMergeResults();
	files =loadFilesAsStructArray('.','_v01.mat');
	etax = [files.eta_x];
	etay = [files.eta_y];
	time = arrayfun(@(x) datenum(x.timestamp),files);
	figure(1);
	plot(etax,etay,'x');
	hold on;
	fplot(@(x) x,[min(etax) max(etax)],'--');
	xlabel('\eta_x');
	ylabel('\eta_y');
	hold off;
	figure(2);
	plot(time,etax,'x');
	hold on;
	plot(time,etay,'ro');
	xlabel('time');
	ylabel('viscosity');
	hold off;
end
