function [path,eta_x,eta_y] = loadFiles(path,patern)

[f,bytes,path]=  dirr(path,patern,'name');
var file;
eta_x = [];
eta_y = [];
for i= 1:length(path)
        disp(strcat('loading ',char(path(i))))
        file = load(char(path(i)));
    if ( isfield(file,'eta_x') && isfield(file,'eta_y') )
        eta_x(i)=file.eta_x;
        eta_y(i)=file.eta_y;
    else
        disp('eta_x does not exist, skip')
    end
    clear file
end
%cellfun(@(x) isfield(x,'eta_x'),file)