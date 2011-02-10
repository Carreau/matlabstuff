function [files] = loadFilesAsStructArray(path,patern)

[f,bytes,path]=  dirr(path,patern,'name');
fileList = arrayfun(@(x) load(char(x)),path,'UniformOutput',false);
files = cellfun(@(x) x,fileList);

%cellfun(@(x) isfield(x,'eta_x'),file)
%arrayfun(@(x) load(char(x)),p)


