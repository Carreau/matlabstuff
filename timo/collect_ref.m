files=dir(['E:\Science\data\matthias\Dataset\*_fit_res.mat']);
    
fit_res_fin=[];
a=[]
cp=[]
for i=1:length(files)
    load(['E:\Science\data\matthias\Dataset\',files(i).name]);
    a=[a,([fit_res.alpha])];
    cp=[cp,([fit_res.cp])];
    if i==1
        fit_res_fin=fit_res;
    else
        fit_res_fin=horzcat(fit_res_fin, fit_res);
    end


end
l(:,1)=cp;
l(:,2)=a;
ls=sortrows(l,1);
di=find(diff(ls(:,1))~=0);