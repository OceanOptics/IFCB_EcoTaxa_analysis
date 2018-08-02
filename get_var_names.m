function objects = get_var_names(matfilename)

% This function can be updated to return any additional desired variables

load(matfilename)

names = fieldnames(keepvar);

obj_status_ind = regexp(names,'\w*object_annotation_status\w*');
obj_cat_ind = regexp(names,'\w*object_annotation_category\w*');
obj_cat_parent_ind = regexp(names,'\w*object_annotation_parent_category\w*');
obj_cat_hier_ind = regexp(names,'\w*object_annotation_hierarchy\w*');
obj_area_pix_ind = regexp(names,'\w*object_summed_area\w*');
obj_vol_pix_ind = regexp(names,'\w*object_summed_biovolume\w*');
obj_ecc_ind = regexp(names,'\w*object_eccentricity\w*');

for ii = 1:numel(names)
    
    if isempty(obj_status_ind{ii}) == 0
        obj_status = getfield(keepvar,names{ii});
        
    elseif isempty(obj_cat_ind{ii}) == 0
        obj_cat = getfield(keepvar,names{ii});
        
    elseif isempty(obj_cat_parent_ind{ii}) == 0
        obj_cat_par = getfield(keepvar,names{ii});
        
    elseif isempty(obj_cat_hier_ind{ii}) == 0
        obj_cat_hier = getfield(keepvar,names{ii});
        
    elseif isempty(obj_area_pix_ind{ii}) == 0
        obj_area_pix = getfield(keepvar,names{ii});
        
    elseif isempty(obj_vol_pix_ind{ii}) == 0
        obj_vol_pix = getfield(keepvar,names{ii});
        
    elseif isempty(obj_ecc_ind{ii}) == 0
        obj_ecc = getfield(keepvar,names{ii});
    end
    
end

% Tara Oceans Polar Circle IFCB data have a slightly different naming
% convenation for summed are and perimeter; check for that here:
if exist('obj_area_pix') == 0
    obj_area_pix_ind = regexp(names,'\w*object_summedarea\w*');
    obj_vol_pix_ind = regexp(names,'\w*object_summedbiovolume\w*');
end

for ii = 1:numel(names)
    
    if isempty(obj_area_pix_ind{ii}) == 0
        obj_area_pix = getfield(keepvar,names{ii});
        
    elseif isempty(obj_vol_pix_ind{ii}) == 0
        obj_vol_pix = getfield(keepvar,names{ii});
    end
    
end

% Need to convert pixels to um
% 3.4 pixels per um (calibration for our IFCB still needed)

% area = cross-sectional area of largest contiguous blob in ROI (pixels^2)
% summed area includes all blobs
obj_area = obj_area_pix./(3.4^2); % in um^2

% volume = volume estimate for the largest blob, assuming cross-sections in the
% third dimension are locally circular (solid of revolution for simple
% shapes; distance map based for complex shapes) Moberg and Sosik 2012 (pixels^3)
% summed volume includes all blobs
obj_vol = obj_vol_pix./(3.4^3); % in um^3

objects.status = obj_status;
objects.cat = obj_cat;
objects.cat_par = obj_cat_par;
objects.cat_hier = obj_cat_hier;
objects.sumarea = obj_area;
objects.sumvol = obj_vol;
objects.ecc = obj_ecc;