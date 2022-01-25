function in_region = is_region_subregion_of(specified_parent_region, subregions)
%%IS_REGION_SUBREGION_OF 
%
%   IN_REGION = is_region_subregion_of(SPECIFIED_REGION, SUBREGIONS)
%   takes a cell array of strings in SUBREGIONS, and looks to see if they
%   are in the region given by SPECIFIED_REGION.
%
%   There is a manually created mapping file `region_mapping.csv` which
%   links subregion strings (e.g. 'VISp5') with a parent region ('VISp').
%   So for example:
%
%       in_region = is_region_subregion_of('VISp', {'CA1', 'VISp2/3', 'VISp5', 'CA1', 'VISp2/3'})
%
%   would return [false, true, true, false, true] in `in_region`.
%
%   Note that the `region_mapping.csv` determines the output of this
%   function.
%
%   Note also that this may not act as you expect, for example:
%
%       in_region = is_region_subregion_of('VISp2/3', {'CA1', 'VISp2/3', 'VISp5', 'CA1', 'VISp2/3'})
%
%   would return [false, false, false, false, false] (unless the
%   `region_mapping.csv` has been edited)
%
%   TODO: use the actual ontology from which the regions were derived
%   instead of a separate csv

mapping = readtable('region_mapping.csv');

in_region = false(length(subregions), 1);

for ii = 1 : length(subregions)
    
    % find index of this subregion in the mapping
    idx = find(strcmp(mapping.subregion, subregions{ii}), 1);
    
    if isempty(idx)
        continue
    end
    
    parent_region = mapping.parent_region{idx};
    
    in_region(ii) = strcmp(parent_region, specified_parent_region);
end
    