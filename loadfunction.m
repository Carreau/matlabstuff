% just a dummy loadfunction which look like doing the same as load('smth')
% but which cn be used in structfun or cellfun.

function [structhandle] = loadfunction(path)
	path;
	structhandle = load(path);
end
	
