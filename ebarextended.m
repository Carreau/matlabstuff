function [skey meanvalue stddevi] = ebarextended(files,averagingField,averagingKeyField)
%  files = loadFilesAsCellArray('.','v01.mat');
  key = m_extractfield(files,averagingKeyField);
  skey  = unique(key);
  clear m r c
  meanvalue = [];
  stddevi = [];
  for k = skey
    [meanvalue(end+1) stddevi(end+1)] =  getMAndSForRate(files,k,averagingField,averagingKeyField);
  end
end

function  [m s] = getMAndSForRate(files,keyvalue,avgfield,keyfield)
  ff = structFilterdByFieldValue(files,keyfield,@(x) x == keyvalue);
  c=cellfun(@(x) getField(x,avgfield),ff);
  s=std(c);
  m=mean(c);
end
