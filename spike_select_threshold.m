function spike_select_threshold(sessions,varargin)

for iSess = 1:length(sessions)
    
    if nargin > 1
        channels = varargin(1);
    else
        channels = lfp_findmda(sessions{iSess});
    end
    
    nChannels = size(channels,2);
    
    thresholds = zeros(1,nChannels);
    outputf = sessions{iSess};
    cd(outputf)
    fullfilename = [outputf,'\matlabData.mat'];
    
    if exist(fullfilename,'file')
        answer = questdlg('You are about to overwrite an existing file. Are you sure you want to continue?', 'Warning!!!', 'No');
        switch lower(answer)
            case 'yes'
                warning('Overwriting previously saved file with new thresholds')
            case 'no'
                display('Skipping session')
                continue
            otherwise
                break
        end
        
    end
    
    
    %Read in timestamps for first 10 seconds of session
    FID = fopen(fullfile(outputf,'timestamps.mda'),'r');
    fseek(FID,0,-1);
    timeStamps = fread(FID,ThSecs*FS,'uint64');
    fclose(FID);
    
    timeStamps=reshape(timeStamps,1,size(timeStamps,1)*size(timeStamps,2));
    
    for iChan = 1:nChannels;
        %Read in samples from current chunk
        datafile = fullfile(sessions{iSess}, [channels{1,iChan}]);
        FID = fopen(datafile,'r');
        fseek(FID,0,-1);
        data = fread(FID,ThSecs*FS,'float32');
        fclose(FID);
        %             data=Nlx2MatCSC(datafile,[0 0 0 0 1],0,4,[timeStamps(chunks(iChunk)), timeStamps(chunks(iChunk+1)-1)]);
        data=reshape(data,1,size(data,1)*size(data,2));
        
        
        plot(timeStamps,data);
        xlim([timeStamps(1) timeStamps(end)])
        
        grid on;
        stringa=sprintf('Select threshold for channel %d',iChan);
        title(stringa,'FontSize',30);
        stringa=sprintf('Insert threshold for channel %d: ',iChan);
        options.Resize='on';
        options.WindowStyle='normal';
        thresholds(iChan)=str2double(cell2mat(inputdlg(stringa,'Threshold selection',1,{'600'},options)));
        close;
        
    end
    
    %Save data
    
    save(fullfilename,'thresholds','-v7.3')
end

end