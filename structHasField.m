function [matchingArray] = structHasField(cellArray,fieldname)
matchingArray = cellfun(@(x) isfield(x,fieldname),cellArray);

