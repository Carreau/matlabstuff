function [matchingArray] = structHasField(arrayArray,fieldname)
matchingArray = arrayfun(@(x) isfield(x,fieldname),arrayArray);

