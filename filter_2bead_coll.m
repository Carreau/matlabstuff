function save_data=filter_2bead_coll(f_path)

%I will give it the preprocesses save_data structure
%FIrst we will filter the data, by a binning, the time resolution should be
%1ms in the end. This will reduce the data
 
 if (nargin == 0)
     [file,spath]=uigetfile('E:\Science\data\collection\elasticity');
    f_path=[spath,filesep,file];
 end
 
% %load the file
 load(f_path)
for i=1:length(save_data)
    %get the effective sampling rate
    es=save_data(i).parameters.Effective_Sampling_Rate.value/2;
    bin_length=round(es/1000); %this ensures a time resolution of about 1ms
    out_data(i).name=save_data(i).name;
    out_data(i).parameters=save_data(i).parameters;
    out_data(i).parameters.Effective_Sampling_Rate.value=es/bin_length;
    out_data(i).d=lin_bin(save_data(i).d,bin_length);
    f_in=save_data(i).f;
    t_in=save_data(i).trap_pos;
    clear f_out t_out;
    for k=1:2
        for l=1:2
            f_out(k,l,:)=lin_bin(f_in(k,l,:),bin_length);
            t_out(k,l,:)=lin_bin(t_in(k,l,:),bin_length);
        end
    end
    out_data(i).f=f_out;
    out_data(i).trap_pos=t_out;        
    
end
save_data=out_data;
%and now we save the dataset with the prefix 'filtered'
save([spath,filesep,'filtered_',file],'save_data');

function out=lin_bin(x,b_size)

%first I ensure that the x is an even multiple of bin size, if not I force
%it by cutting the first values
rest=mod(numel(x),b_size);
if rest~=0
    x(1:rest)=[];
end
%now we reshape and get the mean
x_i=reshape(x,b_size,numel(x)/b_size);
out=mean(x_i,1);
    
