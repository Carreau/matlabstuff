function [files] =  liveViscosity()
	recursiveMergeResults();
	files =loadFilesAsCellArray('.','_v01.mat');
	files = structFilterdByFieldValue(files,'eta_x',@(x) x < 0.10);
	etax = cellfun(@(x) x.eta_x, files);
	etay = cellfun(@(x) x.eta_y, files);
	time = cellfun(@(x) datenum(x.timestamp),files);
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
