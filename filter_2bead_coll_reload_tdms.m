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
% save_data=out_data;
for i=1:length(save_data)
    %now load the TDMS file
    n=i;
    [rdata]=convertTDMS(0,save_data(n).name);
    handles.rdata=rdata;
    handles.dfile=save_data(n).name
    handles=preprocess_data(handles);
    
    a_data.name=save_data(n).name;
    a_data.f=handles.f; plot(squeeze(handles.f(1,1,:)))
    a_data.d=handles.d;
    a_data.parameters=handles.parameters;
    a_data.trap_pos=handles.data(:,1:3,:);

  %  save_data(i)=a_data;
    
    
    %get the effective sampling rate
    es=a_data.parameters.Effective_Sampling_Rate.value/2;
    bin_length=round(es/1000); %this ensures a time resolution of about 1ms
    out_data(i).name=a_data.name;
    out_data(i).parameters=a_data.parameters;
    out_data(i).parameters.Effective_Sampling_Rate.value=es/bin_length;
    out_data(i).d=lin_bin(a_data.d,bin_length);
    f_in=a_data.f;
    t_in=a_data.trap_pos;
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
    

function handles=preprocess_data(handles)
% handles    structure with handles and user data (see GUIDATA)

%here we will preprocess the data. Get absolut positions, and get forces.

rdata=handles.rdata;
%now we extract the data from the TDMS files.
parameters=rdata.Data.Parameters.Root;
for i=1:15
    data1(i,:)=rdata.Data.MeasuredData(3+i).Data;
    data2(i,:)=rdata.Data.MeasuredData(18+i).Data;
end

%now we correct for offsets. This assumes that both beads are traps with
%no force at the beginning!!!!
% This also fixes issues where the first trapwas moved. This ensures that
% from now on always the second trap is moved!
a=find(diff(data2(1,:)+data2(2,:))~=0);
if numel(a)==0
    data_=data2; data2=data1;data1=data_;
    a=find(diff(data2(1,:)+data2(2,:))~=0);
end
data1(4,:)=data1(4,:)-mean(data1(4,1:a(1)));
data1(6,:)=data1(6,:)-mean(data1(6,1:a(1)));
data2(4,:)=data2(4,:)-mean(data2(4,1:a(1)));
data2(6,:)=data2(6,:)-mean(data2(6,1:a(1)));

%additionally, this is the right moment to filter the data, to get rid of
%the annoying noise from the bad digital cable
%ok maybe later, first I just do a simple bin filter:


cal(1)=parameters.AOD_center_X.value;
cal(2)=parameters.AOD_center_Y.value;
cal(3)=parameters.AOD_factor_X.value;
cal(4)=parameters.AOD_factor_Y.value;
cal(5)=parameters.AOD_to_microm_x.value;
cal(6)=parameters.AOD_to_microm_y.value;
xy_slopes(1,1)=parameters.x_slope_bead1.value;
xy_slopes(1,2)=parameters.y_slope_bead1.value;
xy_slopes(2,1)=parameters.x_slope_bead2.value;
xy_slopes(2,2)=parameters.y_slope_bead2.value;

kappa(1,1)=parameters.x_kappa_bead1.value;
kappa(1,2)=parameters.y_kappa_bead1.value;
kappa(2,1)=parameters.x_kappa_bead2.value;
kappa(2,2)=parameters.y_kappa_bead2.value;

%this corrects for the case that the position was saved in digital AOD
%units
if mean(data1(1,:))>1
    data1(1,:)=dig2nor(data1(1,:));
    data1(2,:)=dig2nor(data1(2,:));
    data2(1,:)=dig2nor(data2(1,:));
    data2(2,:)=dig2nor(data2(2,:));
end
    
%recalc the absolute position and put in 12, and 13
data1(12,:)=data1(1,:)/cal(5)-1/xy_slopes(1,1)*data1(4,:)./data1(8,:);
data1(13,:)=data1(2,:)/cal(6)-1/xy_slopes(1,2)*data1(6,:)./data1(8,:);

data2(12,:)=data2(1,:)/cal(5)-1/xy_slopes(2,1)*data2(4,:)./data2(8,:);
data2(13,:)=data2(2,:)/cal(6)-1/xy_slopes(2,2)*data2(6,:)./data2(8,:);

%now we recalc the forces as we want the offset substraction to also act on
%the forces
data1(14,:)=kappa(1,1)/xy_slopes(1,1)*data1(4,:)./data1(8,:)*1E-6;
data1(15,:)=kappa(1,2)/xy_slopes(1,2)*data1(6,:)./data1(8,:)*1E-6;
data2(14,:)=kappa(2,1)/xy_slopes(2,1)*data2(4,:)./data2(8,:)*1E-6;
data2(15,:)=kappa(2,2)/xy_slopes(2,2)*data2(6,:)./data2(8,:)*1E-6;

data(1,:,:)=data1;
data(2,:,:)=data2;

d=sqrt((data1(12,:)-data2(12,:)).^2+(data1(13,:)-data2(13,:)).^2);
%now we change the coordinate system to see the forces parallel and
%perpendicular of the two beads. These forces will be saved in the forces
%array
alp=cart2pol(data1(1,1)-data2(1,1),data1(2,1)-data2(2,1))+pi;
[ft,fr]=cart2pol(data1(14,:),data1(15,:));
[f(1,1,:),f(1,2,:)]=pol2cart(ft+alp,fr);
[ft,fr]=cart2pol(data2(14,:),data2(15,:));
[f(2,1,:),f(2,2,:)]=pol2cart(ft+alp,fr);


handles.d=d;
handles.f=f;
handles.cal=cal;
handles.xy_slopes=xy_slopes;
handles.data=data;
handles.parameters=parameters;

%----------------------This function goes from digital back to normaizer
%AOD units
    function n=dig2nor(d)
        n=(((d/256-2^23)*500/2^23)-75)/15;
