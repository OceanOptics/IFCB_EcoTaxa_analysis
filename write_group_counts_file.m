function [C,IC] = write_group_counts_file(file_loc,matfilename,objects)

obj_val = [];
n = 1;
for ii = 1:length(objects.status)
    if strcmp(objects.status(ii,:), "validated") == 1
        obj_val(n) = ii;
        n = n+1;
    end
end

if isempty(obj_val) == 0
    cell_cat = cellstr(objects.cat);
    [C,~,IC] = unique(cell_cat(obj_val)); % C is groups, IC is the indicies of each group such that A = C(IC)
    counts = accumarray(IC,1);
else 
    C = 'dummy';
    IC = 'dummy2';
end
% Print a table with object counts for each category
ind = strfind(matfilename,'/');
savefilename = matfilename(ind(end)+1:end-4);

fid = fopen([file_loc,'group_counts_text_files/',savefilename,'.txt'],'w');
if isempty(obj_val) == 0
    for ii = 1:length(C)
        fprintf(fid,'%-30s  %-6d\n',C{ii},counts(ii));
    end
end
fprintf(fid,'\n%-30s  %-6d\n','TOTAL VALIDATED',length(obj_val));
fprintf(fid,'\n%-30s  %-6d\n','TOTAL IMAGES',length(objects.status));
fclose(fid);

    