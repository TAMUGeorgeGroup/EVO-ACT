%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% author: Vicky Fan
% The following code is for EMT.
% All the values shown below are arbitury examples. Please feel free to change them by needs.  
% Texas A&M University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if num_tumor_cell_real == EMT_threshould
    new_EMT_cell_ID_set = dividable_tumor_idx;
    EMT_id_history = [EMT_id_history,new_EMT_cell_ID_set];
    D = 100;
    EMT_occurrance_time = n;
elseif num_tumor_cell_real > EMT_threshould
    new_EMT_cell_ID_set = setdiff(dividable_tumor_idx,EMT_id_history);
    EMT_id_history = [EMT_id_history,new_EMT_cell_ID_set];
end

if num_tumor_cell_real >= EMT_threshould
    
    if ~isempty(new_EMT_cell_ID_set)
        tumor_div_rate_percell(new_EMT_cell_ID_set) = tumor_div_rate_EMT;
        tumor_mov_rate_percell(new_EMT_cell_ID_set) = tumor_mov_rate_EMT;
        tumor_div_rate_averaged = sum(tumor_div_rate_percell)/length(tumor_div_rate_percell);
        tumor_mov_rate_averaged = sum(tumor_mov_rate_percell)/length(tumor_mov_rate_percell);
        dividable_cell_clone_id_set = nan(1,length(new_EMT_cell_ID_set));
        for j = 1:length(dividable_cell_clone_id_set)
            dividable_cell_clone_id_set_prep = find(cellfun(@(x) any(x == new_EMT_cell_ID_set(j)), clone_group_ID),1);
            if ~isempty(dividable_cell_clone_id_set_prep)
                dividable_cell_clone_id_set(j) = dividable_cell_clone_id_set_prep;
                dividable_cell_antigen_set = clonegroup_antigen{dividable_cell_clone_id_set(j)};
                if length(dividable_cell_antigen_set) >1
                    if length(dividable_cell_antigen_set) <=15
                        num_antigen_lost_perdividablecell = poissrnd(length(dividable_cell_antigen_set));
                        while num_antigen_lost_perdividablecell >= length(dividable_cell_antigen_set)
                            num_antigen_lost_perdividablecell = poissrnd(length(dividable_cell_antigen_set));
                        end
                    else
                        num_antigen_lost_perdividablecell = poissrnd(15);
                        while num_antigen_lost_perdividablecell >= length(dividable_cell_antigen_set) || num_antigen_lost_perdividablecell ==0
                            num_antigen_lost_perdividablecell = poissrnd(15);
                        end
                    end
                    antigen_to_lost_idx_dividable = randperm(numel(dividable_cell_antigen_set),num_antigen_lost_perdividablecell);
                    dividable_cell_antigen_afterEMT = dividable_cell_antigen_set;
                    dividable_cell_antigen_afterEMT(antigen_to_lost_idx_dividable) =[];
                    clone_group_ID{dividable_cell_clone_id_set(j)} = setdiff(clone_group_ID{dividable_cell_clone_id_set(j)}, new_EMT_cell_ID_set(j)); %Remove the EMT cell from its original clone.
                    clone_group_ID = [clone_group_ID, {new_EMT_cell_ID_set(j)}];
                    clone_group_ID_full{dividable_cell_clone_id_set(j)} = setdiff(clone_group_ID_full{dividable_cell_clone_id_set(j)}, new_EMT_cell_ID_set(j)); %Remove the EMT cell from its original clone.
                    clone_group_ID_full = [clone_group_ID_full, {new_EMT_cell_ID_set(j)}];
                    if isempty(clone_group_ID{dividable_cell_clone_id_set(j)})
                        clonegroup_antigen{dividable_cell_clone_id_set(j)} = nan;
                        averaged_tumor_div_rate_perclone(dividable_cell_clone_id_set(j)) = 0;
                    end
                    clonegroup_antigen = [clonegroup_antigen, {dividable_cell_antigen_afterEMT}];
                    dividable_cell_div_rate_afterEMT = tumor_div_rate_percell(new_EMT_cell_ID_set(j)) + (tumor_div_rate_percell(new_EMT_cell_ID_set(j))/300)*length(dividable_cell_antigen_afterEMT);
                    tumor_div_rate_percell(new_EMT_cell_ID_set(j)) = dividable_cell_div_rate_afterEMT;
                    averaged_tumor_div_rate_perclone = [averaged_tumor_div_rate_perclone, dividable_cell_div_rate_afterEMT];
                end
                parent_daughter_clone(end+1,:) = [dividable_cell_clone_id_set(j), length(clone_group_ID)];
            end
        end
    end
end