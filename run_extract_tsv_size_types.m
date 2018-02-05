function bins_size = run_extract_tsv_size_types(tsv,dilution,write_files,sizebins,diameter)

% Code to process IFCB image data file(s) exported from EcoTaxa
%
%  INPUTS:
%
%   tsv - can either be the name of a single .tsv file from EcoTaxa, or the
%         path of a folder containing multiple .tsv files
%
%   dilution - either a scaler (for one .tsv file) or a vector, set equal to the dilution factor in %; 0 if no dilution 
%
%   write_files - set to 0 to prevent saving a .txt file of group counts 
%               for each sample and a .mat file of analysis results;
%               default is for files to be written and stored in folders as
%               listed below
%
%   sizebins - Can either be a structure defining: number of bins, and min & max 
%               of size bins, which results in calculation of the individual bin
%               widths based on a log scale; OR an array of bin edges
%               (should be length of number of bins you want + 1).
%
%       As a structure of bins number, min, and max:
%
%       sizebins.numbins - the number of size bins
%       sizebins.Dmin - the lower bound of the smallest size bin (um)
%       sizebins.Dmax - the upper bound of the biggest size bin (um)
%
%   diameter - set to 0 for volume-based diameter or 1 for area-based
%               diameter
%
%  OUTPUTS:
%
%   bins_size - a structure containing analysis results. If multiple .tsv
%           files are provided, the data from each file are stored in a structure named
%           bins_file, and bins_size is a structure of structures.
%           For example: biovolume data from the first bin of the first .tsv file is accessed
%           by: bins_size(1).bins_file(1).vol_by_bin
%
%       <struct>.bininds - the indices of all particles found within each bin
%       <struct>.area_by_bin - area of each particle found within a given bin
%       <struct>.vol_by_bin - volume of each particle found within a given bin
%       <struct>.ecc_by_bin - eccentricity of each particle found within a given bin
%       <struct>.category - the EcoTaxa annotation catgory **BOTH PREDICTED AND VALIDATED**
%       <struct>.units - units of area and volume data
%       <struct>.tot_types - total number of images of each type within a given bin
%       <struct>.area_by_type - area of all images of a given type within a given bin
%       <struct>.vol_by_type - volume of all images of a given type within a given bin
%       <struct>.filename - name of origninal .tsv file
%       <struct>.bin_edges - values of bin edges for a given bin (either provided or calculated (um))
%
%  Created folders:
%
%   In the parent folder to the tsv file(s) location, three folders are
%   created to store output files: 
%
%   reduced_mat_files/
%   group_counts_text_files/
%   analysis_mat_files/ 
%
%       If write_files = 0 during function call, the 
%       group_counts_text_files/ and analysis_mat_files/ folders will not
%       be created (reduced_mat_files/ directory is always created).
%
%   Library of required functions:
%   
%   get_var_names.m
%   read_tsv_file.m
%   write_group_counts_file.m
%   bins_by_size.m
%   group_by_type.m
%
% EXAMPLE FUNCTION CALLS
% 
%   EXAMPLE 1: for a folder containing multiple .tsv files, no dilution files (0), yes
%   to have text and results files written (1), a structure of bin size
%   information provided, and volume-based diameter:
%
%       >> tsv = '/Users/lee/Desktop/Leeka/tara/Tara_oceans_polar_circle/IFCB_analysis/ECOTAXA/export_123_20171227_1705_168/';
%       >> sizebins.numbins = 10;
%       >> sizebins.Dmin = 2;
%       >> sizebins.Dmax = 200;
%       >> bins_size = run_extract_tsv_size_types(tsv,0,1,sizebins,0)
%
%   EXAMPLE 2: for a single .tsv file, no dilution (0), no files will be
%   written (0), a vector of bin edges provided, and volume-based diameter:
%
%       >> tsv = 'AT34_INLINE_D20160603T062836_IFCB107.tsv'
%       >> bins_size = run_extract_tsv_size_types(tsv,0,0,[2 5 10 20 100 200],0)
%
%
% A. Chase, University of Maine, Sept 2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%addpath '/Users/lee/Desktop/Leeka/tara/Tara_oceans_polar_circle/IFCB_analysis/ECOTAXA/codes/'

if isdir(tsv) == 0
    file_loc = [];
    files = tsv;
else
    file_loc = tsv;
    files = dir([file_loc,'*.tsv']);
end

% Prepare folders to store output files
if exist([file_loc,'reduced_mat_files'],'dir') ~= 7; mkdir([file_loc,'reduced_mat_files']);end
if write_files ~= 0
    if exist([file_loc,'group_counts_text_files'],'dir') ~= 7; mkdir([file_loc,'group_counts_text_files']);end
    if exist([file_loc,'../analysis_mat_files'],'dir') ~= 7; mkdir([file_loc,'../analysis_mat_files']);end
end

% Loop through all sample files
if isdir(tsv) == 0
    loopnum = 1;
else
    loopnum = length(files);
end

for ifile = 1:loopnum  %For all samples in the "export" folder from EcoTaxa use: 1:loopnum
    
    if isdir(tsv) == 0
        tsvfile = files;
        filename = files;
    else
        tsvfile = files(ifile);  
        filename = tsvfile.name;
    end
    
    disp('***********************************************************')
    disp(['Currently evaluating: ',filename])
    
    % Check to see if the .tsv file from EcoTaxa has been read and reduced yet to
    % removed un-needed variables. If it has, do not repeat (takes time).
    if exist([file_loc,'reduced_mat_files/',filename(1:end-4),'.mat'],'file') == 0
        matfilename = read_tsv_file(tsvfile);
    else
        if isdir(tsv) == 1
            matfilename = [tsvfile.folder,'/reduced_mat_files/',filename(1:end-4),'.mat'];
        else
            matfilename = [file_loc,'reduced_mat_files/',filename(1:end-4),'.mat'];
        end
    end
      
    % Load .mat file of IFCB data and get the relevant variable names
    objects = get_var_names(matfilename);
    
    ids_val = false(length(objects.status),1);
    ids_tot = length(objects.status);
    for ii = 1:ids_tot
        ids_val(ii) = strcmp(objects.status(ii,:), "validated");
    end

    % Print the percent validated to the screen 
    disp([num2str(sum(ids_val)),' out of ',num2str(ids_tot),' images validated (',num2str(round(sum(ids_val)/ids_tot*100,1)),'%)'])
    
    % Find the number of validated objects in each unique category name and print a text file with that info    
    if write_files ~= 0
        [~,~] = write_group_counts_file(file_loc,matfilename,objects); 
    end
       
    % If the sample is from a dilution, multiply the numbers to compare with
    % other per liter data
    if dilution > 0       
        dil_amt = 100./dilution;
        
        area_mat = repmat(objects.sumarea,1,dil_amt);
        objects.sumarea = reshape(area_mat,length(area_mat(1,:))*length(area_mat(:,1)),1);
        
        vol_mat = repmat(objects.sumvol,1,dil_amt);
        objects.sumvol = reshape(vol_mat,length(vol_mat(1,:))*length(vol_mat(:,1)),1);       
    end
      
    % Define bin edges if min, max, and number of bins is provided use
    % sizebins structure
    if isstruct(sizebins) == 1      
        X = (sizebins.Dmax/sizebins.Dmin)^(1/sizebins.numbins);
        bins_right_edge = X.^(1:sizebins.numbins)*sizebins.Dmin;
        
        bin_edges = [sizebins.Dmin,bins_right_edge];       
    else 
        bin_edges = sizebins;
    end
        
    % Size distribution analysis of all particles
    bins = bins_by_size(objects,bin_edges,diameter); 
    
    % Type and group analysis (adds three new fields to the "bins"
    % structure defined just above)
    bins = group_by_type(objects,bin_edges,bins);
    
    % Add the original .tsv filename and bins edges to the structure of binned data
    [bins.filename] = deal(filename);
    for ii = 1:length(bin_edges)-1
        bins(ii).bin_edges = [bin_edges(ii),bin_edges(ii+1)];
    end
    
    % If more than one .tsv file is being analyzed, create a structure of structures.
    if loopnum > 1 
        bins_size(ifile).bins_file = bins;
    else 
        bins_size = bins;
    end      
    
end
    
% Write a .mat file to store the data 
if write_files ~= 0
    getslash = strfind(file_loc,'/');
    save([file_loc,'../analysis_mat_files',file_loc(getslash(end-1):end-1),'_PROC.mat'],'bins_size');
end

    