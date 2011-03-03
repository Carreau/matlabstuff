function Avg = linbin(data,binsize)
    %%
    data=squeeze(data);
    l= size(data,2);
    if (l == 1) && (size(data,1) ~= 1);
        data= data';
    end
    l= size(data,2);
    nRemove = rem(l,binsize);  %# Find the number of points to remov
    data = data(1:end-nRemove);        %# Trim points to make an even multiple of 300
    l= size(data,2);
    %%
    Avg = mean(reshape(data,binsize,l/binsize));

