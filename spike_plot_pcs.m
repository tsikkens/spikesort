function spike_plot_pcs(axhandle,pc,clusters,clus2plot)

axes(axhandle)
cla(axhandle)
hold(axhandle,'on')
grid(axhandle,'on')
title(axhandle,'Projection on the first 2 principal components')

cols_pc=['.k';'.r';'.b';'.m';'.c';'.g';'.r';'.b';'.k';'.m';'.c';'.g';'.r';'.b';'.k';'.m';'.c';'.g';'.r';'.b';'.k';'.m';'.c';'.g';'.k';'.r';'.b';'.m';'.c';'.g';'.r';'.b';'.k';'.m';'.c';'.g';'.r';'.b';'.k';'.m';'.c';'.g';'.r';'.b';'.k';'.m';'.c';'.g'];



for iClus = clus2plot
    
    index = find(clusters == iClus);
    plot(axhandle,pc(index,1),pc(index,2),cols_pc(iClus+1,:));
    
end