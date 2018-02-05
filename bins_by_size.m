function bins = bins_by_size(objects,bin_edges,diameter)

if diameter == 0
% Calculate biovolume-based diameter
    D = ((6./pi)*objects.sumvol).^(1/3); % in um
else
% Calcualte area-based diameter
    D = sqrt(objects.sumarea*4./pi);  % in um
end

%Find the indices of the data within each size bin and populate the
% 'bins' structure (backwards, to reduce memory use)
for ii = length(bin_edges)-1:-1:1
    
    bins(ii).bininds = find(D > bin_edges(ii) & D <= bin_edges(ii+1)); % Indicies of all objects within each bin 
    bins(ii).area_by_bin = objects.sumarea(bins(ii).bininds); % Area of each particle found within a given bin
    bins(ii).vol_by_bin = objects.sumvol(bins(ii).bininds); % Same for volume
    bins(ii).ecc_by_bin = objects.ecc(bins(ii).bininds); % Eccentricity of each particle found within a given bin
    bins(ii).category = objects.cat(bins(ii).bininds,:); % Category of images (**BOTH VALIDATED AND PREDICTED**)
    bins(ii).units = ['area in um^2;',' volume in um^3'];
    
end