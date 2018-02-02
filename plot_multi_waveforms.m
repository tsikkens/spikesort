function plot_multi_waveforms(axhandle,spikes,spikechannel,clusters,clusnums)
cla(axhandle)
hold(axhandle,'on')
nchans = size(spikes,1);
nsamp = size(spikes,3);
scale = max(max(abs(mean(squeeze(spikes(spikechannel,:,:))))))*1.5;
cluscount =0;
nclus = numel(clusnums);


for iclus = 1:nclus;
    index = find(clusters == clusnums(iclus));
    for ichan = 1:nchans
        if ichan == spikechannel
            plot(axhandle,(1:nsamp)+cluscount*nsamp,mean(squeeze(spikes(ichan,index,:)))+scale*(ichan-1),'r')
        else
            plot(axhandle,(1:nsamp)+cluscount*nsamp,mean(squeeze(spikes(ichan,index,:)))+scale*(ichan-1),'k')
        end
    end
    plot(axhandle,[nsamp*cluscount nsamp*cluscount],ylim(axhandle),'k--')
    cluscount = cluscount+1;
end

xlim(axhandle,[0 cluscount*nsamp])