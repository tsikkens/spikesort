function spike_plot_waveforms(axhandle, spikes,clusters,clus2plot)

cols_sp=['-k';'-r';'-b';'-m';'-c';'-g';'-r';'-b';'-k';'-m';'-c';'-g';'-r';'-b';'-k';'-m';'-c';'-g';'-r';'-b';'-k';'-m';'-c';'-g';'-k';'-r';'-b';'-m';'-c';'-g';'-r';'-b';'-k';'-m';'-c';'-g';'-r';'-b';'-k';'-m';'-c';'-g';'-r';'-b';'-k';'-m';'-c';'-g']; % Replace with rand(1,3)?
axes(axhandle)
cla(axhandle)
hold(axhandle,'on')


for iClus = clus2plot
    
index = find(clusters == iClus);

errorbar(mean(squeeze(spikes(index,:))),std(squeeze(spikes(index,:)),0,1),cols_sp(iClus+1,:),'LineWidth',2);

end
xlim([1 64])
legend(axhandle,cellfun(@(x) {num2str(x)}, SplitVec(clus2plot)))