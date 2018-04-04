function spike_plot_isi(axhandle,spikes_ts,clusters,clus2plot)

axes(axhandle)
cla(axhandle)
max_t=50; %ms
edges=0:1:max_t; %Bins (in ms) for plotting ISI distributions
cols_pc=['.k';'.r';'.b';'.m';'.c';'.g';'.r';'.b';'.k';'.m';'.c';'.g';'.r';'.b';'.k';'.m';'.c';'.g';'.r';'.b';'.k';'.m';'.c';'.g';'.k';'.r';'.b';'.m';'.c';'.g';'.r';'.b';'.k';'.m';'.c';'.g';'.r';'.b';'.k';'.m';'.c';'.g';'.r';'.b';'.k';'.m';'.c';'.g'];


%ISI

for iClus = clus2plot
    
    
    
    ts = spikes_ts(clusters == iClus);
    
    isi=(ts(2:end)-ts(1:end-1))/1e3;
    n=histc(isi,edges);
    n=n./sum(n);
    
    bar(edges,n,cols_pc(iClus+1,2));
    xlabel('Time (ms)');
    xlim([edges(1) edges(end)]);
    title('ISI distribution');
    
    
    
end