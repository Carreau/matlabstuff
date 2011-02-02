function [files] =  liveViscosity()
	recursiveMergeResults();
	files =loadFilesAsStructArray('.','_v01.mat');
	etax = arrayfun(@(x) x.eta_x, files);
	etay = arrayfun(@(x) x.eta_y, files);
	time = arrayfun(@(x) datenum(x.timestamp),files);
	figure(1);
	plot(etax,etay,'+');
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
