function [filteredCellArray ] = structWithField(cellArray,fieldName)
ba =  structHasField(cellArray,fieldName);
filteredCellArray = cellArray(ba);
