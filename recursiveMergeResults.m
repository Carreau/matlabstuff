function recursiveMergeResults(patern)
    [f,bytes,p]=  dirr('.','\lts.mat\>','name');
    cellfun(@(x) mergeResultFiles(x),p);
end
