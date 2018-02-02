%% options

%TS Last edit 23-01-2018

datafolder = 'F:\test\PVCre116\Sort'; %data source
% outputfolder = 'E:\Passive_MMN\';%folder to save matlabData.mat file to (creates new subfolder per session)
%
sessions = {datafolder};

chooseThresholds = 0;
ThSecs = 10; %number of seconds to do threshold selection on

nump=64; %Number of samples per waveform
prePoints = 22; %number of samples before spike peak

FS = 32000;


nump=64; %Number of samples per waveform
prePoints = 22; %number of samples before spike peak


for iSess = 1:length(sessions);
    channels = lfp_findmda(sessions{iSess});
    nChannels = size(channels,2);
    outputf = sessions{iSess};
    %     mkdir(outputf)
    cd(outputf)
    fullfilename = [outputf,'\matlabData.mat'];
    MatObj = matfile(fullfilename,'Writable',true);
    
    for ispikeChan = 2:nChannels
        spikes_ind = cell2mat(MatObj.spikes_ind(1,ispikeChan));
        wav = zeros(nChannels,numel(spikes_ind),nump);
        for iChan = 1:nChannels
            
            %Read in samples from current channel
            datafile = fullfile(sessions{iSess}, [channels{1,iChan}]);
            FID = fopen(datafile,'r');
            data = fread(FID,inf,'float32');
            fclose(FID);
            data=reshape(data,1,size(data,1)*size(data,2));
            
            nspikes = size(spikes_ind,2);
            wav_ind = repmat([-prePoints:nump-prePoints-1],nspikes,1);
            sp_ind = repmat([spikes_ind'],1,64);
            index = wav_ind+sp_ind;
            
            wav(iChan,:,:) = reshape(data(index),nspikes,nump);
        end
        
        
        MatObj.spikes_waveforms(1,ispikeChan) = {wav};
    end
    
end