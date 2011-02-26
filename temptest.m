function [s]=temptest(f)
    for i=1:length(f)
        s(i)=twoBeadApprochExperiment(f(i));
    end