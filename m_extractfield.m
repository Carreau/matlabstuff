function fieldstruct =  extractfield(files,field)
	fieldstruct = arrayfun(@(x) getfield(x,field), files);
end
