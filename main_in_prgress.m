%TS Last edit 28-03-2018

foldername = 'E:\Canon_MMN_passiveAwake\1.52 Rec 1\Sort\2018-02-20_11-04-13';
filename = fullfile(foldername,'matlabData.mat');
matObj=matfile(filename,'Writable',true);

% [channels,channelsToAnalyze] = nlx_findncs(foldername);

channelsToAnalyze = 1:size(matObj.spikes_ts,2);
nump=size(cell2mat(matObj.spikes_waveforms(1,1)),3); %Number of samples per waveform

trainingSetSize= 10000;%Inf;
testSetSize = 0.1;
calc_index=1; %Calculate cluster quality indexes? (1=yes, 0=no)

ltol=0.01; %Threshold for the convergence of the gaussian fitting step (normal value 0.01)

cth=0; %Ratio of points to be considered artefacts, set to 0 not to exclude any spike (normal value 0.025 (outermost 5% is excluded in fsmem) or lower)

pcc=2;
th = 0.4;

%Options for the first clustering
opt_smem.Kmin = 1;
opt_smem.Kini = 2; %Used to be 5
opt_smem.Kmax = 10; %Used to be 20
opt_smem.maxite_fullem = 1000; %Used to be 500
opt_smem.maxite_partialem = 1000; %Used to be 500
opt_smem.epsi_fullem = 1e-7; %Used to be 1e-7
opt_smem.epsi_partialem = 1e-7; %Used to be 1e-7
opt_smem.maxcands_split = 12;
opt_smem.maxcands_merge = 12;
opt_smem.fail_exit = 2;

%%
for iChan = channelsToAnalyze
    
    spikes=cell2mat(matObj.spikes_waveforms(1,iChan));
%     spikes= squeeze(spikes(iChan,:,:));
%     
    ts=cell2mat(matObj.spikes_ts(1,iChan));
%     redo_clus = 1;%     
%     while redo_clus
        L=-Inf;
        
        %Prepare test set, to be kept unchanged across multiple runs of the smem algorith         
       
        while L<=0
        testSet=randperm(size(spikes,1));
        testSet=testSet(1:round(size(spikes,1)*testSetSize));             
            
            [clusters,k,pc,L,D,BIC,cont] = SortSpikes_fsmem(spikes,trainingSetSize,pcc,cth,th,ltol,calc_index,0,testSet,opt_smem);
            
            %clusters: vector with the index of the cluster each spikes belongs to (0: artefact cluster)
            %k: number of clusters (+ artefact clusters if present)
            %pc: principal components
            %L: log-likelihood value (the higher it is, the better clustering is)
            %D: Davies-Bouldin's index, useful if clusters were computed on a subset of all spikes (log-likelihood was not computed on the whole data set)
            %Notes:
            %D is not computed (I didn't have the code available, but it can be easily added)
            %BIC: the lower, the better
            %cont: This is the sum of false positives and false negatives per cluster, it must be lower than 0.05, or a cluster must be discarded
        end
%     end

    sp.clusters = clusters;
    sp.nClus = k;
    sp.pc = pc;
    sp.L = L;
    sp.BIC = BIC;
    sp.cont = cont;

    matObj.sp(1,iChan) = {sp};
end