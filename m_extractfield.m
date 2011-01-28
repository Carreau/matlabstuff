function fieldstruct =  extractfield(files,field)
	fieldstruct = cellfun(@(x) getfield(x,field), files);
end
