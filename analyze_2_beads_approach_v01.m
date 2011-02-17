function varargout = analyze_2_beads_approach_v01(varargin)
% ANALYZE_2_BEADS_APPROACH_V0 M-file for analyze_2_beads_approach_v0.fig
%      ANALYZE_2_BEADS_APPROACH_V0, by itself, creates a new ANALYZE_2_BEADS_APPROACH_V0 or raises the existing
%      singleton*.
%
%      H = ANALYZE_2_BEADS_APPROACH_V0 returns the handle to a new ANALYZE_2_BEADS_APPROACH_V0 or the handle to
%      the existing singleton*.
%
%      ANALYZE_2_BEADS_APPROACH_V0('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ANALYZE_2_BEADS_APPROACH_V0.M with the given input arguments.
%
%      ANALYZE_2_BEADS_APPROACH_V0('Property','Value',...) creates a new ANALYZE_2_BEADS_APPROACH_V0 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before analyze_2_beads_approach_v0_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to analyze_2_beads_approach_v0_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help analyze_2_beads_approach_v0

% Last Modified by GUIDE v2.5 01-Nov-2010 01:33:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @analyze_2_beads_approach_v0_OpeningFcn, ...
                   'gui_OutputFcn',  @analyze_2_beads_approach_v0_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before analyze_2_beads_approach_v0 is made visible.
function analyze_2_beads_approach_v0_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to analyze_2_beads_approach_v0 (see VARARGIN)

% Choose default command line output for analyze_2_beads_approach_v0
handles.output = hObject;

%DEfile a starting directory fpr the open file
handles.ddir='E:\Science\data\elasticity';

%make sure we have the tools
set(hObject,'toolbar','figure');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes analyze_2_beads_approach_v0 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = analyze_2_beads_approach_v0_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function load_Callback(hObject, eventdata, handles)
% hObject    handle to load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ddir=handles.ddir;
[dfile,ddir]=uigetfile('*.tdms','Please specify the TDMS file to analyze',ddir);
handles.actual_set=1;
file_list.name=[ddir,filesep,dfile];
handles.file_list=file_list
handles.ddir=ddir;
handles.dfile=dfile;
%now make sure that the save data array is empty
handles.save_data=[];
load_dataset(hObject, eventdata, handles)



% --------------------------------------------------------------------
function handles=preprocess_data(handles)
% handles    structure with handles and user data (see GUIDATA)

%here we will preprocess the data. Get absolut positions, and get forces.

%r is the bead radius, unfortunately, I do not save that one from the
%preogram
r=4.5e-6/2
%this is the force over which I try to align all the datasets
f_align=10e-12

rdata=handles.rdata;
%now we extract the data from the TDMS files.
parameters=rdata.Data.Parameters.Root;


%% filter data
% due to control bit we need to remove 1 point from time to time
% let's find the middle value and find on which side ther is the
% most points
clear rd milieu want s;
rd = rdata.Data.MeasuredData(4).Data;
milieu  = (min(rd)+max(rd))/2;
want = rd < milieu;
s= sum(want);
if(s < length(rd)/10)
   want=not(want);
end
clear s milieu;
disp(sprintf('in the filtered data, %d datapoints have been removed', sum(not(want)) ));

%% try to change the rdata structure to remove unwanted columns
rd2=rdata
rd=rdata;
arrayfun(@(x) disp(length(x.Data)),rd2.Data.MeasuredData)
s = sum(not(want))
for i = 1:33
    rdata.Data.MeasuredData(i).Data = rdata.Data.MeasuredData(i).Data(want);
    rdata.Data.MeasuredData(i).Total_Samples = rdata.Data.MeasuredData(i).Total_Samples-sum(not(want));
end
disp('après filtrage');
arrayfun(@(x) length(x.Data),rd2.Data.MeasuredData)
%rdata.Data.MeasuredData


%%

disp('fin filtragee');
for i=1:15
    data1(i,:)=rdata.Data.MeasuredData(3+i).Data;
    data2(i,:)=rdata.Data.MeasuredData(18+i).Data;
end
%%%


%now we correct for offsets. This assumes that both beads are traps with
%no force at the beginning!!!!
% This also fixes issues where the first trap help was moved. This ensures that
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

%now we will already filter the data to have it less heavy, and also apply
%the calculation of the values of interest
    a_data.name=handles.file_list(handles.actual_set).name;
    a_data.f=handles.f; plot(squeeze(handles.f(1,1,:)))
    a_data.d=handles.d;
    a_data.parameters=handles.parameters;
    a_data.trap_pos=handles.data(:,1:3,:);
  
    %get the effective sampling rate
    es=a_data.parameters.Effective_Sampling_Rate.value/2;
    bin_length=round(es/1000); %this ensures a time resolution of about 1ms
    out_data.name=a_data.name;
    out_data.parameters=a_data.parameters;
    out_data.parameters.Effective_Sampling_Rate.value=es/bin_length;
    out_data.d=lin_bin(a_data.d,bin_length);
    f_in=a_data.f;
    t_in=a_data.trap_pos;
    clear f_out t_out;
    for k=1:2
        for l=1:2
            f_out(k,l,:)=lin_bin(f_in(k,l,:),bin_length);
            t_out(k,l,:)=lin_bin(t_in(k,l,:),bin_length);
        end
    end
    out_data.f=f_out;
    out_data.trap_pos=t_out;  
    

  %and now we calculate the numbers we would like to extract from the data
      %first I calculate the dissipated energy
    out_data.diss_energy=1e-6*sum(diff(out_data.d)'.*squeeze(out_data.f(1,1,2:end)));
    set(handles.text_diss,'String',['E Diss=',num2str(out_data.diss_energy,2),' J'])
    
    %now get teh points where the approach stated and ended, and where the
    %retraction started
    tr=out_data.trap_pos;
    tr_m=sqrt(squeeze(tr(2,1,:)).^2+squeeze(tr(2,2,:)).^2);
    tr_d_0=find(diff(tr_m)==0);
    tr_d_no0=find(diff(tr_m)~=0);
    
    out_data.appr_start=tr_d_no0(1);
    out_data.retr_start=tr_d_0(end);
    tr_d_0_cut=tr_d_0;
    tr_d_0_cut(find(tr_d_0_cut<=tr_d_no0(1)))=[];
    out_data.appr_stop=tr_d_0_cut(1);
    
    %so now I need to estimate the 'touching point'. I will try to do this
    %by getting the point where the forces are bigger than f_align
    [fm,pm]=min(abs(f_align-squeeze(out_data.f(1,1,1:out_data.appr_stop))));
    out_data.touch_point=pm;
   
    %I will also try to get the touching point for the retraction phase 'touching point'. I will try to do this
    %by getting the point where the forces are bigger than f_align
    [fm,pm_r]=min(abs(f_align-squeeze(out_data.f(1,1,out_data.retr_start:end))));
    pm_r=pm_r+out_data.retr_start;
    out_data.touch_point_retr=pm_r;
    
    
    %This point is now used to calculate the indentation
    out_data.indent=out_data.d(pm)-out_data.d(out_data.retr_start);
    set(handles.text_indent,'String',['Indentation=',num2str(out_data.indent,2),' µm'])
    
    %Now we can also use the 'touching point as reference for the elasticy
    %estimate as a function of the indentation. The formular is:
    %F=4/3*E/(1-(nu)^2).*d.^(3/2)*sqrt(r)+f0;
    f_ind=squeeze(out_data.f(1,1,pm:out_data.appr_stop))';
    d_ind=out_data.d(pm:out_data.appr_stop);
    d_ind=-(d_ind-d_ind(1))*1e-6;
    out_data.E=(f_ind-f_align)*3*(1-0.5^2)./(4*d_ind.^(3/2)*sqrt(r));
    out_data.E_final=mean(out_data.E(round(numel(out_data.E)/2:end)));

    set(handles.text_E,'String',['Youngs Mod=',num2str(out_data.E_final,2),' Pa'])   
    
    
    %Next is the exponential decay first extract the regime from the data,
    %then fit an exponential of the form d(t)=alpha*exp(-t/tau)+d_e
    
    d_t=out_data.d(out_data.appr_stop:out_data.retr_start);
    t_t=[0:length(d_t)-1]*1/out_data.parameters.Effective_Sampling_Rate.value;
    
    start_point = [d_t(1)-d_t(end),d_t(end),t_t(end)];
    %options=optimset('Display','iter');
    %estimates = fminsearch(@exp_dec, start_point,options,t_t,d_t);
    estimates = fminsearch(@exp_dec, start_point,[],t_t,d_t);
    [sse,exp_fit]=exp_dec(estimates,t_t,d_t);
    estimates;
    out_data.alpha=estimates(1);
    out_data.tau=estimates(3);
    out_data.d0=estimates(2);
    out_data.rel_decay=estimates(1)/estimates(2);
    
    set(handles.text_tau,'String',['Rel Time=',num2str(out_data.tau,2),' s'])       
    
    %finally I will get the right regimes for the overlay replotting
    out_data.f_app_overlay=squeeze(out_data.f(1,1,1:out_data.appr_stop));
    out_data.f_retr_overlay=squeeze(out_data.f(1,1,out_data.retr_start:end));
    out_data.d_app_overlay=out_data.d(1:out_data.appr_stop)-out_data.d(pm);
    out_data.d_retr_overlay=out_data.d(out_data.retr_start:end)-out_data.d(pm_r);

    
    %finally we put it all in the handles structure to be saved later if
    %necessary
    handles.out_data=out_data;
%guidata(hObject, handles);
% --------------------------------------------------------------------
function update_figures(hObject,handles)
% handles    structure with handles and user data (see GUIDATA)

%here we will update the plots, including all, data, selection, fits
%First we check which range was chosen. If no range was chosen so far, we take all the data
pos=cell2mat(get(findobj(handles.plot_ft,'Tag','separator'),'XData'));
if length(pos)==0
    pos=[.01 .01; ones(1,2)*(length(squeeze(handles.f(1,1,:)))-1)*2/(handles.parameters.Effective_Sampling_Rate.value)];
end

%make sure that the first slider is for a smaller value than the second
if diff(pos(:,1))<0
    pos=flipud(pos);
end
%check if any of the values is out of bounds
pos(1,1)=inrange(pos(1,1),0,(length(squeeze(handles.f(1,1,:)))-1)*2/(handles.parameters.Effective_Sampling_Rate.value));
pos(2,1)=inrange(pos(2,1),0,(length(squeeze(handles.f(1,1,:)))-1)*2/(handles.parameters.Effective_Sampling_Rate.value));



axes(handles.plot_fd)
sel1=round(pos(1,1)/2*(handles.parameters.Effective_Sampling_Rate.value))+1;
sel2=round(pos(2,1)/2*(handles.parameters.Effective_Sampling_Rate.value));
d_sel=handles.d(sel1:sel2);
f_sel=squeeze(handles.f(1,1,sel1:sel2));
plot(d_sel,f_sel);
title('Force Distance Plot')
xlabel('Distance between beads [µm]')
ylabel('Force in N')

axes(handles.plot_ft)
ha=gca;
%the position of the actual selector lines
plot([0:length(squeeze(handles.f(1,1,:)))-1]*2/(handles.parameters.Effective_Sampling_Rate.value),squeeze(handles.f(1,1,:)));
hold on;
plot([0:length(squeeze(handles.f(2,1,:)))-1]*2/(handles.parameters.Effective_Sampling_Rate.value),squeeze(handles.f(2,1,:)),'k')
hold off;

title('Force Evolution Plot')
xlabel('Time [s]')
ylabel('Force in N')
%now we check if we already have the lines, if not, we create them
if length(cell2mat(get(findobj(handles.plot_ft,'Tag','separator'),'XData')))==0
    selectgui(hObject,handles,ha,2,pos)
end

axes(handles.plot_ft_sel)
ha=gca;
%the position of the actual selector lines
plot([0:length(d_sel)-1]*2/(handles.parameters.Effective_Sampling_Rate.value),f_sel);
title('Force Evolution Plot')
xlabel('Time [s]')
ylabel('Force in N')

%update the handles structure
handles.d_sel=d_sel';
handles.f_sel=f_sel;

guidata(hObject, handles);

%  ---- This is the function from Tobias, that creates the two lines
function selectgui(hObject,handles,ha, nos, varargin)
% SELECTGUI     Creates multiple dragable vertical separator lines in a given axis
%
% SELECTGUI(HA, NOS) creates NOS dragable vertical separator
% lines in axis with handle HA.
%
% SELECTGUI(HA, NOS, STARTPOS) does the same as above but
% additionally, the vector STARTPOS defines the Startposition of the
% separator lines. The length of the vector STARTPOS has to be equal to
% the number of separators (NOS).
%
% If STARTPOS is not provided, the separator lines are linearly spaced.
%
% Usage:
%   1) do your plottings in any axis with hande HA
%   2) call selectgui (separator lines will be created)
%   3) at any time, you can get the current x-positions of the 
%       separator lines by calling
%
%               positions = cell2mat(get(findobj(HA,'Tag','separator'),'XData'));
%
%   4) If you don't need the separators anymore just call
%
%               delete(findobj(HA,'Tag','separator'));
%
%       or simply close the window.
%   
%   5) Have fun !
%
% Example:
%
%       ha = gca;                                   % handle to the current axis
%       data_x = rand(100,1);                       % plot some random data
%       data_y = rand(100,1);
%       plot(ha, data_x, data_y, '+');
%
%       logpos = logspace(min(data_x), max(data_x),7) ./ 10;    % prepare a logarithmically scaled vector
%       selectgui(ha, 5, logpos(2:6));              % call selectgoi
%
%  The positions of the separators can be received at any time by calling
%
%       positions = cell2mat(get(findobj(ha,'Tag','separator'),'XData'));
%
%  The separators can be deleted by calling
%
%       delete(findobj(ha,'Tag','separator'));
%
% % %  % %  % %  % %  % %  % %  % %  % %  % %  % %  % %  % %  % %  % %
% % (c) Tobias Kießling, University of Leipzig, Germany, Oct.21.2010 %
% % %  % %  % %  % %  % %  % %  % %  % %  % %  % %  % %  % %  % %  % %


axisState = get(ha, 'NextPlot');    % Store the current state of the 'NextPlot' property of the axis

set(ha, 'NextPlot', 'add');             % the same as 'hold on'
hf = get(ha, 'parent');                 % handle to the figure
set(hf, 'WindowButtonUpFcn', {@vsep_callback_buttonUp,hObject,handles});  % set the WindowButtonUpFcn of the figure
ylim = get(ha,'ylim');                  % y - limits of the axis
xlim = get(ha,'xlim');                  % x - limits of the axis


% % create separatorLines
x = repmat(struct('start',[]),nos,1);
for i = 1 : nos
    if isempty(varargin)                                                        % Did the user specify a startvector ?
        x(i).start = xlim(1) + (i)/(nos+1) * diff(xlim);                        % -> No!  generate equally spaced separator lines
    else
        x(i).start = varargin{1}(i);                                                    % -> Yes! so use the passed starting positions
    end

    hsep = plot(ha, [x(i).start, x(i).start], [ylim(1), ylim(2)],'r');          % plot the separator line
    set(hsep, 'Tag', 'separator');                                                    % the separator line gets the Tag 'separator'
    set(hsep, 'ButtonDownFcn', @vsep_callback_buttonDown)           % ascripe ButtonDownFcn to separator line
end
% Lines created

set(ha, 'NextPlot', axisState);    % restore the intital 'NextPlot' property


%these are the callbacks for bars

    function vsep_callback_buttonDown(hsep, varargin,handles)       % is called when the user clicks on a separator 
                                                                                     % (hsep is the handle to the separator)
        ha = get(hsep, 'parent');       % get the handle to the axis
        hfigure = get(ha, 'parent');    % get the handle to the figure

        setappdata(ha, 'HSelSep', hsep);    % store the handle to the selected separator within the figure
        set(hfigure, 'WindowButtonMotionFcn', @vsep_callback_moving)    % Enable WindowButtonMotionFcn
    


    function vsep_callback_buttonUp(hfigure, varargin,hObject,handles)
        set(hfigure, 'WindowButtonMotionFcn', '')           % Button is up -> Disable WindowButtonMotionFcn
        update_figures(hObject,handles);


function vsep_callback_moving(hfigure, varargin)    % The WindowsButtonMotionFcn:
        ha = get(hfigure, 'CurrentAxes');                           % get the handle to the current axis
        mousepos = get(ha, 'CurrentPoint');                      % get the current position of the mouse pointer
        HSelSep = getappdata(ha, 'HSelSep');                   % get the handle of the current separator line
                                                                                 % (the handle was stored in the figure. See buttonDownFcn)

        set(HSelSep, 'XData', [mousepos(1), mousepos(1)]);  % Move the seperator line to the current mouseposition
    
    function x=inrange(x,low,up)
        if x<low
            x=low;
        elseif x> up
            x=up
        end


% --- Executes on selection change in fit_sel.
function fit_sel_Callback(hObject, eventdata, handles)
% hObject    handle to fit_sel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns fit_sel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from fit_sel

%here we will do the fits for the data, depending on the choise
c=get(handles.fit_sel,'Value');
if c==1
    handles=fit_fd(hObject,handles);
elseif c==2
    handles=fit_ft(handles);
end

% --- Executes during object creation, after setting all properties.
function fit_sel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fit_sel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% -------- Do the force distanc fit
function handles=fit_fd(hObject,handles);

d_sel=handles.d_sel;
f_sel=handles.f_sel;
r=4.5e-6/2;
n=.5;

%prepare the f and d, as the fit function is for the case of first touch
d=(max(d_sel)-d_sel(1:10:end))*1e-6;
f=f_sel(1:10:end)-min(f_sel);



start_point = [1000 0.01];
options=optimset('Display','iter');
estimates = fminsearch(@herz_function, start_point,options,d,f,r,n);

[sse,f_fit]=herz_function(estimates,d,f,r,n);
axes(handles.plot_fd)
hold on
plot(-d*1e6+max(d_sel),f_fit+min(f_sel),'r');
set(handles.text4,'String',['Youngs modulus E [Pa]= ',num2str(estimates(1))]);
hold off
estimates(1)
    
    
% -------- Do the force evolution fit
function handles=fit_ft(handles);

f_sel=handles.f_sel;


%prepare the f and d, as the fit function is for the case of first touch
t=[0:length(f_sel)-1]*2/(handles.parameters.Effective_Sampling_Rate.value);
t=t(1:10:end);
f=f_sel(1:10:end);

start_point = [f(1)-f(end),f(end),t(end)];
options=optimset('Display','iter');
estimates = fminsearch(@exp_dec, start_point,options,t',f);

[sse,exp_fit]=exp_dec(estimates,t',f);
axes(handles.plot_ft_sel)
hold on
plot(t,exp_fit,'r');
set(handles.text5,'String',['Relaxation time tau [s]= ',num2str(estimates(3))]);
hold off
estimates(3)
    
    
% ----------------------------------------------------------------------    
%the function for the fmin search on the force distance
function [sse, FittedCurve] =herz_function(params,d,f,r,nu)
        %global eta2 fc
        E = params(1);
        f0 = params(2);
        FittedCurve=4/3*E/(1-(nu)^2).*d.^(3/2)*sqrt(r)+f0;
        ErrorVector = (FittedCurve - f)./f;
        ErrorVector(isnan(ErrorVector))=[];
        sse = abs(sum(abs(ErrorVector) .^ 2));
  
        
% ----------------------------------------------------------------------    
%the function for the fmin search on the force distance
function [sse, FittedCurve] =exp_dec(params,t,f)
        %global eta2 fc
        a = params(1);
        b = params(2);
        tau= params(3);
        FittedCurve=a*exp(-t./tau)+b;
        ErrorVector = (FittedCurve - f)./f;
        ErrorVector(isnan(ErrorVector))=[];
        sse = abs(sum(abs(ErrorVector) .^ 2))*1e10;        
        
        
%----------------------This function goes from digital back to normaizer
%AOD units
    function n=dig2nor(d)
        n=(((d/256-2^23)*500/2^23)-75)/15;


% --- Executes on button press in take_data.
function take_data_Callback(hObject, eventdata, handles)
% hObject    handle to take_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%now we store_the actuall_data_and go to the next file
% a_data.name=handles.file_list(handles.actual_set).name;
% a_data.f=handles.f;
% a_data.d=handles.d;
% a_data.parameters=handles.parameters;
% a_data.trap_pos=handles.data(:,1:3,:);

a_data=handles.out_data;

if length(handles.save_data)==0
    handles.save_data=a_data;
else
    handles.save_data(end+1)=a_data;
end
set(handles.num_in_save_data,'String',['The current number in the save_data struct is ',num2str(length(handles.save_data))]);

%This checks if the actual set is already the last. If not, it moves one on
n=handles.actual_set;
if n+1<=length(handles.file_list)
    handles.actual_set=n+1;
end
load_dataset(hObject, eventdata, handles)
%guidata(hObject, handles);

% --- Executes on button press in ignore_data.
function ignore_data_Callback(hObject, eventdata, handles)
% hObject    handle to ignore_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
n=handles.actual_set;
if n+1<=length(handles.file_list)
    handles.actual_set=n+1;
end
load_dataset(hObject, eventdata, handles)

% --------------------------------------------------------------------
function load_dir_Callback(hObject, eventdata, handles)
% hObject    handle to load_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ddir=handles.ddir;
ddir=uigetdir(ddir,'Please specify the parent dir');
%now make sure that the save data array is empty
handles.save_data=[];
set(handles.num_in_save_data,'String','The current number in the save_data struct is 0');
file_list=get_subfolders(ddir);
handles.actual_set=1;
handles.file_list=file_list;
handles.ddir=ddir;
guidata(hObject, handles);

%this now loads the next dataset
load_dataset(hObject, eventdata, handles)




%-----------------------------------------------------------
function file_list=get_subfolders(parent_dir)
%disp(parent_dir);
allSubFolders=genpath(parent_dir);
    
    % Scan through them separating them.
remain = allSubFolders;
listOfFolderNames = {};
n=1;
while true % Demo code adapted from the help file.
[singleSubFolder, remain] = strtok(remain, pathsep);
if isempty(singleSubFolder), break; end
disp(sprintf('%s', singleSubFolder));
listOfFolderNames = [listOfFolderNames singleSubFolder];
d_t=dir([listOfFolderNames{end},filesep,'approach_2beads*.tdms']);
    for i=1:length(d_t)
        file_list(n).name=[listOfFolderNames{end},filesep,d_t(i).name];
        n=n+1;
    end
end


function load_dataset(hObject, eventdata, handles)
%now load the TDMS file
n=handles.actual_set;
file_list=handles.file_list;
[rdata]=convertTDMS(0,file_list(n).name);
set(handles.text3,'String',file_list(n).name)

handles.rdata=rdata;
handles.dfile=file_list(n).name;
guidata(hObject, handles);

handles=preprocess_data(handles);
update_figures(hObject,handles);


guidata(hObject, handles);


% --------------------------------------------------------------------
function save_dat_Callback(hObject, eventdata, handles)
% hObject    handle to save_dat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%now we ask the user where we should store the whole stuff
%[sfile,spath]=uigetfile('Where should I store the save_data',handles.ddir)
save_data=handles.save_data;
uisave('save_data');

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
