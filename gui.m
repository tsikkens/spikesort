function varargout = gui(varargin)
% GUI MATLAB code for gui.fig
%      
% See also: GUIDE, GUIDATA, GUIHANDLES

% Last Modified by GUIDE v2.5 04-Apr-2018 16:54:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @gui_OpeningFcn, ...
    'gui_OutputFcn',  @gui_OutputFcn, ...
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


% --- Executes just before gui is made visible.
function gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui (see VARARGIN)

% Choose default command line output for gui
global MAT sp iChan
MAT = [];
sp = [];
iChan = [];

handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIW
% AIT makes gui wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles;


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1
global MAT iChan sp

listStr = get(hObject,'String');
listVal = get(hObject,'Value');


spikes =  cell2mat(MAT.spikes_waveforms(1,iChan));
spikes_ts = cell2mat(MAT.spikes_ts(1,iChan));
clus2plot = unique(sp.clusters);

if listVal ~=1
    clusSelect = str2double(listStr{listVal});
    clus2plot = clusSelect;
end

spike_plot_waveforms(handles.ProbeMapAxes,spikes,sp.clusters,clus2plot)
spike_plot_isi(handles.IsiAxes,spikes_ts,sp.clusters,clus2plot)

nspikes = sum(ismember(sp.clusters,clus2plot));
stringm='False negative/positive indexes:';
for j=clus2plot
    if j == 0
        continue
    end
    stringm=sprintf('%s\n Cluster %d: %d',stringm,j,sp.cont(j));
end
text = {['Channel: ' num2str(iChan)];['nSpikes: ' num2str(nspikes)];stringm};
set(handles.clusterData,'String',text)


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in SplitButton.
function SplitButton_Callback(hObject, eventdata, handles)
% hObject    handle to SplitButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global sp iChan MAT

listStr = get(handles.listbox1,'String');
listVal = get(handles.listbox1,'Value');


if listVal ~= 1
    indx = find(sp.clusters == str2double(listStr{listVal}));
else
    warndlg('Please select a single cluster')
end

sp_new = cluster_split_function(MAT,iChan,indx);
sp.clusters(indx) = sp_new.clusters + max(sp.clusters);
sp.cont = [sp.cont sp_new.cont];

spikes =  cell2mat(MAT.spikes_waveforms(1,iChan));
spikes_ts = cell2mat(MAT.spikes_ts(1,iChan));
clus2plot = unique(sp.clusters);

liststr = {'all'};
for iClus = clus2plot
    liststr = cat(1,liststr,num2str(iClus));
end


set(handles.listbox1,'Value',1)
set(handles.listbox1,'String',liststr)

clus2plot(clus2plot == 0) = [];



spike_plot_waveforms(handles.WaveformAxes,spikes,sp.clusters,clus2plot)
spike_plot_pcs(handles.PcAxes,sp.pc,sp.clusters,clus2plot)
spike_plot_isi(handles.IsiAxes,spikes_ts,sp.clusters,clus2plot)

nspikes = sum(ismember(sp.clusters,clus2plot));
stringm='False negative/positive indexes:';
for j=clus2plot
    stringm=sprintf('%s\n Cluster %d: %d',stringm,j,sp.cont(j));
end
text = {['Channel: ' num2str(iChan)];['nSpikes: ' num2str(nspikes)];stringm};
set(handles.clusterData,'String',text)






% --- Executes on button press in MergeButton.
function MergeButton_Callback(hObject, eventdata, handles)
% hObject    handle to MergeButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global sp iChan MAT

R = inputdlg('Which clusters do you want to merge?','Merge',1,{'[]'});
R = eval(cell2mat(R));

sp.clusters(sp.clusters == R(2)) = R(1);

spikes =  cell2mat(MAT.spikes_waveforms(1,iChan));
spikes_ts = cell2mat(MAT.spikes_ts(1,iChan));
clus2plot = unique(sp.clusters);

liststr = {'all'};
for iClus = clus2plot
    liststr = cat(1,liststr,num2str(iClus));
end


set(handles.listbox1,'Value',1)
set(handles.listbox1,'String',liststr)

clus2plot(clus2plot == 0) = [];



spike_plot_waveforms(handles.WaveformAxes,spikes,sp.clusters,clus2plot)
spike_plot_pcs(handles.PcAxes,sp.pc,sp.clusters,clus2plot)
spike_plot_isi(handles.IsiAxes,spikes_ts,sp.clusters,clus2plot)

nspikes = sum(ismember(sp.clusters,clus2plot));
stringm='False negative/positive indexes:';
for j=clus2plot
    stringm=sprintf('%s\n Cluster %d: %d',stringm,j,sp.cont(j));
end
text = {['Channel: ' num2str(iChan)];['nSpikes: ' num2str(nspikes)];stringm};
set(handles.clusterData,'String',text)





% --- Executes on button press in AssignToNoiseButton.
function AssignToNoiseButton_Callback(hObject, eventdata, handles)
% hObject    handle to AssignToNoiseButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global sp iChan MAT
listStr = get(handles.listbox1,'String');
listVal = get(handles.listbox1,'Value');

if listVal ~= 1
    sp.clusters(sp.clusters == str2double(listStr{listVal})) = 0;
else
    sp.clusters(:) = 0;
end

spikes =  cell2mat(MAT.spikes_waveforms(1,iChan));
spikes_ts = cell2mat(MAT.spikes_ts(1,iChan));
clus2plot = unique(sp.clusters);
liststr = {'all'};
for iClus = clus2plot
    liststr = cat(1,liststr,num2str(iClus));
end


set(handles.listbox1,'Value',1)
set(handles.listbox1,'String',liststr)

clus2plot(clus2plot == 0) = [];

spike_plot_waveforms(handles.WaveformAxes,spikes,sp.clusters,clus2plot)
spike_plot_pcs(handles.PcAxes,sp.pc,sp.clusters,clus2plot)
spike_plot_isi(handles.IsiAxes,spikes_ts,sp.clusters,clus2plot)

nspikes = sum(ismember(sp.clusters,clus2plot));
stringm='False negative/positive indexes:';
for j=clus2plot
    stringm=sprintf('%s\n Cluster %d: %d',stringm,j,sp.cont(j));
end
text = {['Channel: ' num2str(iChan)];['nSpikes: ' num2str(nspikes)];stringm};
set(handles.clusterData,'String',text)





% --- Executes on button press in ConfirmClustersButton.
function ConfirmClustersButton_Callback(hObject, eventdata, handles)
% hObject    handle to ConfirmClustersButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global MAT iChan sp

spikes_clusters = sp.clusters;
clustersToKeep = unique(spikes_clusters);
count = 0;
for iClus = clustersToKeep
    if iClus ==0
        continue
    end
    count = count +1;
    spikes_clusters(spikes_clusters == iClus) = count;
end

clustersToKeep = unique(spikes_clusters);
clustersToKeep(clustersToKeep == 0) = [];
MAT.spikes_clusters(1,iChan) = {spikes_clusters};
MAT.clustersToKeep(1,iChan) = {clustersToKeep};

iChan = iChan +1;
sp = cell2mat(MAT.sp(1,iChan));
spikes =  cell2mat(MAT.spikes_waveforms(1,iChan));
spikes_ts = cell2mat(MAT.spikes_ts(1,iChan));
clus2plot = unique(sp.clusters);
liststr = {'all'};
for iClus = clus2plot
    liststr = cat(1,liststr,num2str(iClus));
end
set(handles.listbox1,'Value',1)
set(handles.listbox1,'String',liststr)
clus2plot(clus2plot == 0) = [];
spike_plot_waveforms(handles.WaveformAxes,spikes,sp.clusters,clus2plot)
spike_plot_pcs(handles.PcAxes,sp.pc,sp.clusters,clus2plot)
spike_plot_isi(handles.IsiAxes,spikes_ts,sp.clusters,clus2plot)

nspikes = sum(ismember(sp.clusters,clus2plot));
stringm='False negative/positive indexes:';
for j=clus2plot
    stringm=sprintf('%s\n Cluster %d: %d',stringm,j,sp.cont(j));
end
text = {['Channel: ' num2str(iChan)];['nSpikes: ' num2str(nspikes)];stringm};
set(handles.clusterData,'String',text)




% --- Executes on button press in loaddatabutton.
function loaddatabutton_Callback(hObject, eventdata, handles)
% hObject    handle to loaddatabutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global MAT iChan sp
[filename, pathname] = uigetfile('*matlabData.mat','Select matlabData file');
ffile = fullfile(pathname,filename);
MAT = matfile(ffile,'Writable',true);
iChan = 1;
sp = cell2mat(MAT.sp(1,iChan));
spikes =  cell2mat(MAT.spikes_waveforms(1,iChan));
spikes_ts = cell2mat(MAT.spikes_ts(1,iChan));
clus2plot = unique(sp.clusters);
liststr = {'all'};
for iClus = clus2plot
    liststr = cat(1,liststr,num2str(iClus));
end


set(handles.listbox1,'Value',1)
set(handles.listbox1,'String',liststr)

clus2plot(clus2plot == 0) = [];
spike_plot_waveforms(handles.WaveformAxes,spikes,sp.clusters,clus2plot)
spike_plot_pcs(handles.PcAxes,sp.pc,sp.clusters,clus2plot)
spike_plot_isi(handles.IsiAxes,spikes_ts,sp.clusters,clus2plot)

nspikes = length(sp.clusters);
stringm='False negative/positive indexes:';
for j=1:length(sp.cont)
    stringm=sprintf('%s\n Cluster %d: %d',stringm,j,sp.cont(j));
end
text = {['Channel: ' num2str(iChan)];['nSpikes: ' num2str(nspikes)];stringm};
set(handles.clusterData,'String',text)


% --- Executes on button press in undoButton.
function undoButton_Callback(hObject, eventdata, handles)
% hObject    handle to undoButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global MAT iChan sp
sp = cell2mat(MAT.sp(1,iChan));
spikes =  cell2mat(MAT.spikes_waveforms(1,iChan));
spikes_ts = cell2mat(MAT.spikes_ts(1,iChan));
clus2plot = unique(sp.clusters);
liststr = {'all'};
for iClus = clus2plot
    liststr = cat(1,liststr,num2str(iClus));
end


set(handles.listbox1,'Value',1)
set(handles.listbox1,'String',liststr)

clus2plot(clus2plot == 0) = [];
spike_plot_waveforms(handles.WaveformAxes,spikes,sp.clusters,clus2plot)
spike_plot_pcs(handles.PcAxes,sp.pc,sp.clusters,clus2plot)
spike_plot_isi(handles.IsiAxes,spikes_ts,sp.clusters,clus2plot)

nspikes = length(sp.clusters);
stringm='False negative/positive indexes:';
for j=1:length(sp.cont)
    stringm=sprintf('%s\n Cluster %d: %d',stringm,j,sp.cont(j));
end
text = {['Channel: ' num2str(iChan)];['nSpikes: ' num2str(nspikes)];stringm};
set(handles.clusterData,'String',text)


% --- Executes during object deletion, before destroying properties.
guide

% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
