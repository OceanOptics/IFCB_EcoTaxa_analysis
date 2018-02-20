function bins = group_by_type(objects,bin_edges,bins)

% Function to add the fields "area_by_type" and "vol_by_type" to the bins
% structure

% This function is currently assigning images to the following types:
% Ciliophora, Diatoms, Dinoflagellates, Dictyochales, and Prymnesiophytes

% Additionally, there are grouped categories for:
% non-living, other, multiple, and unvalidated

% All categories, and their numbers:
rem_cat = [];   %1  to be removed: bubble, bead, artefact, badfocus, part 

non_liv = [];   %2  detritus, feces, plastic, fiber
cili = [];      %3  Ciliates
chlor = [];     %4  Chlorophytes (includes Prasinophytes)
crypt = [];     %5  Cryptophytes
diats = [];     %6  Diatoms
dinos = [];     %7  Dinoflagellates
dict = [];      %8  Dictyota
eug = [];       %9  Euglenophytes
prym = [];      %10 Prymnesiophytes (includes phaeocystis and coccolithophores)
other = [];     %11 other & othertocheck
multiple = [];  %12 
clumps = [];    %13

obj_unval = []; %14 all unvalidated images

obj_numbers = 1:14;

% Future:
% Group together all categories with fewer than 1% of all validated images from the
% current sample
% Add "unidentifiable" category

% Define the non-living categories 
remove_cat = {'\w*bubble\w*','\w*bead\w*','\w*artefact\w*','\w*badfocus\w*','\w*part\w*'};

nonliv_cat = {'\w*detritus\w*','\w*feces\w*','\w*fiber\w*','\w*plastic\w*'};
 
ciliates = {'\w*Ciliophora\w*'};
   
chloros = {'\w*Chlorophytes\w*','\w*Prasinophyceae\w*','\w*tempPrasinophyceae\w*'};

cryptos = {'\w*tempCryptophyceae\w*'};

% Define the diatoms
diatoms = {'\w*pennate\w*','\w*centric\w*','\w*chain\w*','\w*Corethron\w*',...
    '\w*Pseudo-nitzschia\w*','\w*Bacteriastrum\w*','\w*Planktoniella\w*',...
    '\w*Thalassiosira\w*','\w*Coscinodiscus\w*','\w*Membraneis\w*',...
    '\w*Chaetoceros\w*','\w*Guinardia\w*','\w*Rhizosolenia\w*','\w*Navicula\w*',...
    '\w*Ditylum\w*','\w*Bacteriastrum\w*','\w*Cylindrotheca\w*','\w*Eucampia\w*'};

% Define the dinoflagellates
dinoflagellates = {'\w*Dinophyceae\w*','\w*Pyrocystis\w*','\w*Ceratium\w*',...
    '\w*Dinophysis\w*','\w*Oxytoxum\w*','\w*Prorocentrum\w*','\w*Warnowia\w*',...
    '\w*Nematopsides\w*'};

dictyos = {'\w*Dictyo\w*'}; % for example you could add : "t002" here

euglenos = {'\w*Euglen\w*'};

prymnesio = {'\w*Prymnesiophyceae\w*','\w*Phaeocystis\w*','\w*Prymnesiaceae X\w*'};

others = {'\w*other\w*','\w*othertocheck\w*'};%,'\w*t0\w*'};

mult = {'\w*multiple\w*'};

clum = {'\w*clumps\w*'};

% Not used currently:
%     "Tintinnida","nauplii","Mollusca",...
%     "Coccolithales","Phaeocystis","Rhizaria","Foraminifera"
%     (and others)

% *** No further group definition done below this line ***
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find the indicides of validated objects
obj_val = zeros(length(objects.status),1);
n = 1;
for ii = 1:length(objects.status)
    if strcmp(objects.status(ii,:), "validated") == 1 || strcmp(objects.status(ii,:),'"validated"') == 1
        obj_val(n) = ii;
        n = n+1;
    end
end

obj_unval = setdiff(1:1:length(objects.cat),obj_val);

n = 1; m = 1; p = 1; q = 1; r = 1; s = 1; t = 1; u = 1; v = 1; w = 1; x = 1; y = 1; z = 1; 
for ii = 1:length(objects.cat)
    
    % only categorize validated objects
    if sum(obj_val == ii) == 1
        
        for jj = 1:length(remove_cat)     
            if regexp(objects.cat(ii,:), remove_cat{jj}) ~= 0
                rem_cat(n) = ii;
                n = n+1;
            end
        end
        
        for jj = 1:length(nonliv_cat)     
            if regexp(objects.cat(ii,:), nonliv_cat{jj}) ~= 0
                non_liv(m) = ii;
                m = m+1;
            end
        end
        
        for jj = 1:length(ciliates)     
            if regexp(objects.cat(ii,:), ciliates{jj}) ~= 0
                cili(p) = ii;
                p = p+1;
            end
        end
        
        for jj = 1:length(chloros)     
            if regexp(objects.cat(ii,:), chloros{jj}) ~= 0
                chlor(q) = ii;
                q = q+1;
            end
        end
        
        for jj = 1:length(cryptos)
            if regexp(objects.cat(ii,:), cryptos{jj}) ~= 0
                crypt(r) = ii;
                r = r+1;
            end
        end
        for jj = 1:length(diatoms)
            if regexp(objects.cat(ii,:), diatoms{jj}) ~= 0
                diats(s) = ii;
                s = s+1;
            end
        end
        
        for jj = 1:length(dinoflagellates)
            if regexp(objects.cat(ii,:), dinoflagellates{jj}) ~= 0
                dinos(t) = ii;
                t = t+1;
            end
        end
        
        for jj = 1:length(dictyos)
            if regexp(objects.cat(ii,:), dictyos{jj}) ~= 0
                dict(u) = ii;
                u = u+1;
            end
        end
        
        for jj = 1:length(euglenos)
            if regexp(objects.cat(ii,:), euglenos{jj}) ~= 0
                eug(v) = ii;
                v = v+1;
            end
        end
        
        for jj = 1:length(prymnesio)
            if regexp(objects.cat(ii,:), prymnesio{jj}) ~= 0
                prym(w) = ii;
                w = w+1;
            end
        end
        
        for jj = 1:length(others)
            if regexp(objects.cat(ii,:), others{jj}) ~= 0
                other(x) = ii;
                x = x+1;
            end
        end
        
        for jj = 1:length(mult)
            if regexp(objects.cat(ii,:), mult{jj}) ~= 0
                multiple(y) = ii;
                y = y+1;
            end
        end
        
        for jj = 1:length(clum)
            if regexp(objects.cat(ii,:), clum{jj}) ~= 0
                clumps(z) = ii;
                z = z+1;
            end
        end        
    end 
end

% Assign a type number to the index of each object
obj_array = 1:1:length(objects.cat);

obj_array(rem_cat)   = obj_numbers(1);
obj_array(non_liv)   = obj_numbers(2);
obj_array(cili)      = obj_numbers(3);
obj_array(chlor)     = obj_numbers(4);
obj_array(crypt)     = obj_numbers(5);
obj_array(diats)     = obj_numbers(6);
obj_array(dinos)     = obj_numbers(7);
obj_array(dict)      = obj_numbers(8);
obj_array(eug)       = obj_numbers(9);
obj_array(prym)      = obj_numbers(10);
obj_array(other)     = obj_numbers(11); %be sure to change the number of other in the loop below
obj_array(multiple)  = obj_numbers(12);
obj_array(clumps)    = obj_numbers(13);
obj_array(obj_unval) = obj_numbers(14);

% Put any remaining objects that do not fall into one of the groups into
% the other category and print a message to the screen
obj_remain = find(obj_array > max(obj_numbers));
if length(obj_remain) >= 1
    obj_array(obj_remain) = obj_numbers(11);
    disp(['There were ',num2str(length(obj_remain)),' objects moved to "other" after not fitting in a defined category, their types are:']) 
    for ii = 1:length(obj_remain)
        disp(objects.cat(obj_remain(ii),:)) 
    end
end

% Print the number of objects to the screen that will be removed in the
% loop below
rem_inds = find(obj_array==1);
disp(['There will be ',num2str(length(rem_inds)),' objects removed from the analysis'])

objects.typenum = obj_array';

% Calculate the type information for each bin previously determined with
% size analysis
for ii = length(bin_edges)-1:-1:1
        
    bins_types = objects.typenum(bins(ii).bininds); % Indicies of each type in a given bin
    bins_sumarea = objects.sumarea(bins(ii).bininds); % Area of particles in a given bin
    bins_sumvol = objects.sumvol(bins(ii).bininds); % Volume of particles in a given bin 

    keep_inds = find(bins_types ~= 1);
    if isempty(keep_inds) == 1
        keep_inds = zeros(0,1);
    end
    
    bins(ii).tot_types = accumarray(bins_types(keep_inds),1); % Accumulates over all types indexed by 1 through the highest number in bins(ii).types
    
    bins(ii).area_by_type = accumarray(bins_types(keep_inds),bins_sumarea(keep_inds)); % Summed area of all particles of each type within a given bin
    bins(ii).vol_by_type = accumarray(bins_types(keep_inds),bins_sumvol(keep_inds)); % Same for volume
    
    % Populate with zeros if accumarray did not find values for the categories
    % at the end of the list
    if isempty(bins_types(keep_inds)) == 1 % condition for when there are no validated particles in a bin
        length_diff = max(obj_numbers);
    else
        length_diff = abs(max(bins_types(keep_inds)) - max(obj_numbers));
    end
    
    if length_diff > 0
        
        bins(ii).tot_types = [bins(ii).tot_types',zeros(1,length_diff)]';
        bins(ii).area_by_type = [bins(ii).area_by_type',zeros(1,length_diff)]';
        bins(ii).vol_by_type = [bins(ii).vol_by_type',zeros(1,length_diff)]';
        
    end
    
    bins(ii).types = 'non-living = 2, cilliates = 3, chlorophytes = 4, cryptophytes = 5, diatoms = 6, dinoflagellates = 7, dictyota = 8, euglenoids = 9, prymnesiophytes = 10, other = 11, multiple = 12, clumps = 13, unvalidated = 14';

end