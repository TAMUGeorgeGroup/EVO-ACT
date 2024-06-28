%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% author: Vicky Fan
% The following code is for tumor mutation.
% All the values shown below are arbitury examples. Please feel free to change them according to needs.  
% Texas A&M University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if rand <= mutation_rate
    clone_group_ID = [clone_group_ID,{num_tumor_cell_track}]; 
    clone_group_ID_full = [clone_group_ID_full,{num_tumor_cell_track}];
    if length(clone_group_ID) <=10
        tumor_counts_per_clone(length(clone_group_ID),n) = 1;
    end
    coin = rand;
    if coin < 0.5 || length(parent_cell_antigen) == 1 
        number_antigen_to_add = 0;
        while number_antigen_to_add==0 || number_antigen_to_add > (length(antigen_set) - length(parent_cell_antigen))
            number_antigen_to_add = poissrnd(10);
        end
        
        antigen_to_add = [];
        while length(antigen_to_add) < number_antigen_to_add
            new_antigen_set = setdiff(antigen_set,[parent_cell_antigen,antigen_to_add]);
            if isempty(new_antigen_set)
                antigen_to_add = [];
            else
                new_antigen = new_antigen_set(randperm(length(new_antigen_set),1));
                antigen_to_add = [antigen_to_add, new_antigen];
            end
        end
        
        daughter_cell_antigen = [parent_cell_antigen, antigen_to_add]; 
        
        if rand <0.5
            new_cell_div_rate = tumor_div_rate_percell(index_div_cell) - number_antigen_to_add/100;
            denom = 101;
            while new_cell_div_rate<=0
                new_cell_div_rate = tumor_div_rate_percell(index_div_cell) - number_antigen_to_add/denom;
                denom = denom +1;
            end
        else
            new_cell_div_rate = tumor_div_rate_percell(index_div_cell) + number_antigen_to_add/1.5;
        end
    elseif coin >= 0.5 || length(parent_cell_antigen) == length(antigen_set)
        number_antigen_to_lost = 0;
        while number_antigen_to_lost==0 || number_antigen_to_lost >= length(parent_cell_antigen)
            number_antigen_to_lost = poissrnd(10);
        end
        antigen_to_lost_idx = randperm(numel(parent_cell_antigen),number_antigen_to_lost);
        daughter_cell_antigen = parent_cell_antigen;
        daughter_cell_antigen(antigen_to_lost_idx) =[];
        if rand <0.5
            new_cell_div_rate = tumor_div_rate_percell(index_div_cell) - number_antigen_to_lost/100;
            denom = 101;
            while new_cell_div_rate<=0
                new_cell_div_rate = tumor_div_rate_percell(index_div_cell) - number_antigen_to_lost/denom;
                denom = denom +1;
            end
        else
            new_cell_div_rate = tumor_div_rate_percell(index_div_cell) + number_antigen_to_lost/1.5;
        end
    end
    clonegroup_antigen = [clonegroup_antigen, {daughter_cell_antigen}];
    tumor_div_rate_percell(num_tumor_cell_track) = new_cell_div_rate;
    averaged_tumor_div_rate_perclone = [averaged_tumor_div_rate_perclone, new_cell_div_rate];
    parent_daughter_clone(end+1,:) = [parent_cell_clone,length(clone_group_ID)];
    
else 
    for i = 1:length(clone_group_ID)
        if any(clone_group_ID{i} == index_div_cell)
            clone_group_ID{i} = [clone_group_ID{i}, num_tumor_cell_track];
            break
        end
    end
    
    for i = 1:length(clone_group_ID_full)
        if any(clone_group_ID_full{i} == index_div_cell)
            clone_group_ID_full{i} = [clone_group_ID_full{i}, num_tumor_cell_track];
            break
        end
    end
    if parent_cell_clone <=10
        tumor_counts_per_clone(parent_cell_clone,n) = tumor_counts_per_clone(parent_cell_clone,find(~isnan(tumor_counts_per_clone(parent_cell_clone,:)), 1, 'last')) + 1;
    end
end
