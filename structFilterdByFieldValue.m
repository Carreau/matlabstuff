function [filteredCellArray,swf,fieldvalue] = structFilteredByFieldValue(cellarray,fieldname,filteringfunction)
swf = structWithField(cellarray,fieldname);
fieldvalue = cellfun(@(x) getfield(x,fieldname),swf,'UniformOutput',false);
f= filteringfunction;
boolarray =  cellfun(f,fieldvalue);
filteredCellArray = swf(boolarray);
