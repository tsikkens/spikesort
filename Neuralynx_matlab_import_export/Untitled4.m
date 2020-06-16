% determine task order in recording session based on events file

filename = 'E:\mmn\1.57 Rec 3\Events.nev';


% load in events file
[TimeStamps, TTLs, EventStrings, Header] = ...
    Nlx2MatEV( filename, [1 0 1 0 1], 1, 1, []);

% check what the original filename is
[og_fname,~] = regexp(Header(2,:),'....-..-.._..-..-..','match','split');
og_fname = og_fname{1};

% determine start and stop of each recording
startrec = find(strcmp(EventStrings,'Starting Recording'));
stoprec = find(strcmp(EventStrings,'Stopping Recording'));

% determine number of recordings
nrecs = length(startrec);

% for each recording determine start and stop timestamps
start_ts = TimeStamps(startrec);
stop_ts = TimeStamps(stoprec);


%% for each recording determine task type 

task = cell(1,nrecs);

for iRecording = 1:nrecs
   
    tasktype = [];
    values = TTLs(startrec(iRecording):stoprec(iRecording));
    values = values(values ~= 0 & values ~= 31 & values ~= 1234);
    
    
    
    uni = unique(values);
    
    if ismember(36,uni)
        tasktype = 'optotagging';
    end
    
    d = diff(values);
    duni = unique(d);
    
    if ismember(1:8,uni)
        iscontrol = true;
    else 
        iscontrol = false;
    end
    
    if iscontrol
        if d(1:6) ~= 0
            tasktype = 'MS';
        else
            tasktype = 'DO';
        end
    else
       if ismember(1:4,uni)
           tasktype = 'AV_MMN';
       elseif ismember([1, 2],uni)
           tasktype = 'AO_MMN';
       elseif ismember([1, 2],uni)
           tasktype = 'VO_MMN';
       end
    end
    
    if isempty(tasktype)
        tasktype = 'baseline/remove';
    end
    
    task{iRecording} = tasktype;
    
    
end