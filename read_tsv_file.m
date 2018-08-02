function [matfilename] = read_tsv_file(tsvfile)

% Read the data exported from EcoTaxa (.tsv files) and remove the unecessary
% variables to make a more managable .mat file
% A. Chase March 2017
filename = tsvfile.name;
filefold = tsvfile.folder;

disp('*** Reducing .tsv file ***')

allvar = tdfread([filefold,'/',filename]);

names = fieldnames(allvar);
n = 1;
for ii = 1:numel(names)
    if isempty(regexp(names{ii},'\w*wedge\w*|\w*_hog\w*|\w*_ring\w', 'once')) == 0
        remfield(n) = ii;
        n = n+1;
    end
end
keepvar = rmfield(allvar,names(remfield));

matfilename = [filefold,'/reduced_mat_files/',filename(1:end-4),'.mat'];

save(matfilename)