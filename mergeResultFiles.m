% function to merge 2 types of files, 
% calibration_simons10XX_AOM=0_bead=N_results.mat
% and
% calibration_simons10XX_AOM=0_bead=N_psd.mat
% use a patern matching % calibration_simons10XX_AOM=0_bead=N_results.mat ('*results.mat')
% and it will look for the psd corresponding file and merge the 2 in a '..merged_v01.mat' file

function mergeResultsFiles(patern)
	d = dir(patern);
	l= length('results.mat');
	names_orig		= arrayfun(@(x) x.name			,d			,'UniformOutput',false);
	base_names		= cellfun( @(x) (x(1:end-l))	,names_orig	,'UniformOutput',false);
	names_extended	= cellfun( @(x) [x 'psd.mat']	,base_names	,'UniformOutput',false);
	names_newoutput	= cellfun( @(x) [x 'merged_v01.mat'],base_names	,'UniformOutput',false);
	cellfun(@(x,y,z) actualmerge(x,y,z),names_orig,names_extended,names_newoutput);
end


function actualmerge(resultfilename,extendedfilname,outputfilename)
	x = loadfunction(resultfilename);
	y = loadfunction(extendedfilname);
	c = catstruct(x,y,'sorted');
	save(outputfilename,'-struct','c');
end

