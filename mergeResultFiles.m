% function to merge 2 types of files, 
% calibration_simons10XX_AOM=0_bead=N_results.mat
% and
% calibration_simons10XX_AOM=0_bead=N_psd.mat
% use a patern matching % calibration_simons10XX_AOM=0_bead=N_results.mat ('*results.mat')
% and it will look for the psd corresponding file and merge the 2 in a '..merged_v01.mat' file

function mergeResultsFiles(filename)
	l= length('results.mat');
	base_names		= filename(1:end-l);
	names_extended	= [base_names 'psd.mat'];
	names_newoutput	= [base_names 'merged_v01.mat'];
    actualmerge(filename,names_extended,names_newoutput);
%	cellfun(@(x,y,z) actualmerge(x,y,z),filename,names_extended,names_newoutput);
end

function actualmerge(resultfilename,extendedfilname,outputfilename)
    disp(['merging into ' outputfilename ]);
	x = loadfunction(resultfilename);
	y = loadfunction(extendedfilname);
	c = catstruct(x,y,'sorted');
    c.filename=outputfilename;
    c.makedAsDefective = 0;
    c.defectiveReason = '';
    c.markAsDefectiveTimeStamp = 0;
	save(outputfilename,'-struct','c');
end

