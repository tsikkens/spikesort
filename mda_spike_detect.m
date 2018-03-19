%% options

%TS Last edit 23-01-2018

datafolder = 'E:\Canon_MMN_passiveAwake\1.52 Rec 1\Sort\2018-02-20_11-04-13'; %data source
% outputfolder = 'E:\Passive_MMN\';%folder to save matlabData.mat file to (creates new subfolder per session)
%
sessions = {datafolder};

chooseThresholds = 0;
ThSecs = 10; %number of seconds to do threshold selection on

nump=64; %Number of samples per waveform
prePoints = 22; %number of samples before spike peak

FS = 32000;


%% Threshold selection
tic
runtime = zeros(1,length(sessions));

%Select all thresholds
if chooseThresholds
    spike_select_threshold(sessions,ThSecs)
end

%% Main

for iSess = 1:length(sessions)
    tic
    
    channels = lfp_findmda(sessions{iSess});
    nChannels = size(channels,2);
    
    %     outputf =fullfile(outputfolder, sessions{iSess});
    outputf = sessions{iSess};
    %     mkdir(outputf)
    cd(outputf)
    fullfilename = [outputf,'\matlabData.mat'];
    MatObj = matfile(fullfilename,'Writable',true);
    
    try %check if thresholds are predefined
        load(fullfilename,'thresholds');
        thFlag = true;
        if ~exist('thresholds','var')
            thFlag = false;
        end
    catch
        thFlag = false;
    end
    
    for iShank = 1%:size(channels,1)
        
        %Read in timestamps and Fs for entire session
        dataset = fullfile(sessions{iSess}, [channels{iShank,1} '.ncs']);
        
        %         nSamples = size(timeStamps,2);
        %         nChannels = size(channels,2);
        %
        
        
        %Preallocate all variables to be stored in matlabData.mat to avoid
        %extreme file sizes
        %
        %         spikes_ts=cell(1,nChannels);
        %         spikes_waveforms=cell(1,nChannels);
        
        h = waitbar(0,'Running Spike Detection...');
        
        for iChan = 1:nChannels
            
            %Read in samples from current channel
            datafile = fullfile(sessions{iSess}, [channels{1,iChan}]);
            FID = fopen(datafile,'r');
            fseek(FID,0,-1);
            data = fread(FID,inf,'float32');
            fclose(FID);
            data=reshape(data,1,size(data,1)*size(data,2));
            
            % Threshold selection
            if thFlag
                th=thresholds(iChan);
            else
                th = mean(data)+4*std(data);
            end
            
            % Spike extraction
            
            if th>0
                app=find(data>th);
                app1=SplitVec(app,'consecutive');
                %                     app=[];
                app = nan(1,length(app1));
                for j=1:length(app1)
                    [maxvalue, maxind]=max(data(1,app1{j}));
                    app(j)=app1{j}(maxind);
                end
            elseif th<0
                app=find(data<th);
                app1=SplitVec(app,'consecutive');
                %                     app=[];
                app = nan(1,length(app1));
                for j=1:length(app1)
                    [minvalue, minind]=min(data(1,app1{j}));
                    app(j)=app1{j}(minind);
                end
            end
            %Check for spikes too close to the beginning or end of
            %chunk
            if ~isempty(app)
                while app(1)<nump
                    app(1)=[];
                    if isempty(app)
                        break;
                    end
                end
                if ~isempty(app)
                    while (length(data)-app(end))<nump
                        app(end)=[];
                        if isempty(app)
                            break;
                        end
                    end
                end
            end
            
            
            MatObj.spikes_ind(1,iChan) = {app};
            
            waitbar(iChan/nChannels,h);
        end
        
        
    end
    clear data
    
    FID = fopen(fullfile(outputf,'timestamps.mda'),'r');
    timeStamps = fread(FID,inf,'uint64');
    fclose(FID);
    timeStamps=reshape(timeStamps,1,size(timeStamps,1)*size(timeStamps,2));
    
    MatObj.spikes_ts = cellfun(@(x) {timeStamps(x)},MatObj.spikes_ind);
    
    close(h)
    
    runtime(iSess) = toc;
    sprintf('Runtime for session: %s was %d seconds',sessions{iSess},runtime(iSess))
end

