function [] = f_get_topo(UniqueVals_sorted,UniqueFreqs,RF)
load('~/Documents/montage_CSD.mat'); %this is the channel list

    rf_idx = find(UniqueFreqs == RF);
    delta_scc = squeeze(nanmedian(UniqueVals_sorted(rf_idx,:,:),2));
    chans_to_plot = 1:64;
    chans_to_plot(ismember(chans_to_plot,[5,13,19])) =[];
%     ismember(chans_to_plot,[5,13,19]) =[];
    delta_scc(find(delta_scc==0)) =[];
    f_topo(E(chans_to_plot),delta_scc)
caxis([-.8,.8])
end 
