function [s m] = ebar(files)
%  files = loadFilesAsCellArray('.','v01.mat');
  rate = m_extractfield(files,'rate');
  srate  = unique(rate);
  clear m r c
  m = [];
  s= [];
  for r = srate
    [s(end+1) m(end+1)] =  getMAndSForRate(files,r);
  end
end

function  [m s] = getMAndSForRate(files,trate)
  ff = structFilterdByFieldValue(files,'rate',@(x) x == trate);
  c=cellfun(@(x) x.eta_x,ff);
  s=std(c);
  m=mean(c);
end