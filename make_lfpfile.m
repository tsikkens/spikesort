clear all;
close all;

%% Options
datafolder = 'F:\Passive_MMN'; %data source
outputfolder = 'F:\test';%folder to save matlabData.mat file to (creates new subfolder per session)

sessions = {'PVCre116\Sort'};

channels = {'CSC1', 'CSC2' , 'CSC3', 'CSC4',  ...
    'CSC5', 'CSC6', 'CSC7', 'CSC8'};%,    ...
%     'CSC9', 'CSC10', 'CSC11', 'CSC12', ...
%     'CSC13', 'CSC14', 'CSC15', 'CSC16',...
%     'CSC17','CSC18','CSC19','CSC20',  ...
%     'CSC21','CSC22','CSC23','CSC24',   ...
%     'CSC25','CSC26','CSC27','CSC28',  ...
%     'CSC29','CSC30','CSC31','CSC32'};%, ...  % probe 1
%     'CSC33', 'CSC34', 'CSC35', 'CSC36', ...
%     'CSC37', 'CSC38', 'CSC39', 'CSC40', ...
%     'CSC41', 'CSC42', 'CSC43', 'CSC44', ...
%     'CSC45', 'CSC46', 'CSC47', 'CSC48', ...
%     'CSC49', 'CSC50', 'CSC51', 'CSC52', ...
%     'CSC53', 'CSC54', 'CSC55', 'CSC56', ...
%     'CSC57', 'CSC58', 'CSC59', 'CSC60', ...
%     'CSC61', 'CSC62', 'CSC63', 'CSC64'}; % probe 2

nChunkSamp = 1024000; %Divide data into chunks of 'nChunkSamp' samples to avoid memory issues. Should be a multiple of 512 to avoid problems



doReref = 1;

doFiltering = 1;
freqLow = 700; %frequency for high-pass filtering
freqHigh = 9000; %frequency for low-pass filtering
filtOrd = 4; %order of butterworth filter used


%% Main
tic
for iSess = 1:length(sessions)
    
    outputf =fullfile(outputfolder, sessions{iSess});
    if ~exist(outputf,'dir')
        mkdir(outputf)
    end
    %Read in timestamps and Fs for entire session
    dataset = fullfile(datafolder,sessions{iSess}, [channels{1,1} '.ncs']);
    [timeStamps, Fs] = Nlx2MatCSC(dataset,[1 0 1 0 0],0,1);
    
    Fs=Fs(1);
    
    app=(0:1E6/Fs:511E6/Fs)';
    app=repmat(app,1,length(timeStamps));
    timeStamps=repmat(timeStamps,512,1);
    timeStamps=timeStamps+app;
    timeStamps=reshape(timeStamps,1,size(timeStamps,1)*size(timeStamps,2));
    
    nSamples = size(timeStamps,2);
    nChannels = size(channels,2);
    
    %Divide session data into chunks of 'nChunkSamp' samples to avoid memory issues
    chunks = [1:nChunkSamp:nSamples,nSamples+1];
    
    FID = fopen(fullfile(outputf,'timestamps.mda'),'w');
    csc_ts=reshape(timeStamps,size(timeStamps,2),1);
    fwrite(FID,csc_ts,'uint64');
    fclose(FID);
    
    for iShank = 1%:size(channels,1)
        h = waitbar(0,'Please wait...');
        
        
        for iChan = 1:nChannels
            stringa=sprintf('matlabData%d.mda',iChan);
            fullfilename{iChan} = fullfile(outputf,stringa);
            FF(iChan)=fopen(fullfilename{iChan},'w');
        end
        
        for iChunk = 1:length(chunks)-1
            
            
            waitbar(iChunk/(length(chunks)-1),h);
            
            if iChunk < length(chunks)-1
                csc_data = zeros(nChannels,nChunkSamp);
            else
                csc_data = [];
            end
            
            for iChan = 1:nChannels
                
                %Read in samples from current chunk
                datafile = fullfile(datafolder,sessions{iSess}, [channels{iShank,iChan} '.ncs']);
                data=Nlx2MatCSC(datafile,[0 0 0 0 1],0,4,[timeStamps(chunks(iChunk)), timeStamps(chunks(iChunk+1)-1)]);
                data=reshape(data,1,size(data,1)*size(data,2));
                if doFiltering
                    [B,A]=butter(filtOrd,[freqLow freqHigh]/(Fs/2));
                    csc_data(iChan,:)=filtfilt(B,A,data);
                else
                    csc_data(iChan,:)=data;
                end
            end
            if doReref
                data_mean = repmat(mean(csc_data),nChannels,1);
                csc_data = csc_data-data_mean;
                clear data_mean;
            end
            
            
            for iChan = 1:nChannels
                csc_data_app=single(reshape(csc_data(iChan,:),size(csc_data,2),1));
                fwrite(FF(iChan),csc_data_app,'float32');
            end
            
        end
        close(h);
        
    end
end

for iChan = 1:nChannels
    fclose(FF(iChan));
end
runtime = toc;
