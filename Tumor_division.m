%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% author: Vicky Fan
% The following code is for tumor division.
% All the values shown below are arbitury examples. Please feel free to change them by needs.  
% Texas A&M University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_steps = 2000000;
Time = nan(1,time_steps);
Time(1) = 0;
real_time = nan(1,time_steps);
real_time(1) = 0;
for n = 2:time_steps
    tumor_inv_waiting = rand./(tumor_inv_rate * num_tumor_pheri);
    tumor_div_waiting = rand./(tumor_div_rate_averaged * num_tumor_cell_real);
    tumor_mov_waiting = rand./(tumor_mov_rate_averaged * num_tumor_cell_real);
    T_move_waiting = rand./(T_move_rates * num_T);
    
    t = [tumor_inv_waiting tumor_div_waiting tumor_mov_waiting T_move_waiting];
    t_new = min(t);
    
    Time(n) = t_new;
    real_time(n) = real_time(n-1) + Time(n);
    index = find(t==t_new);
    num_T_real_pertime(n) = num_T_real_pertime(n-1);
    num_T_recognizable(n) = num_T_recognizable(n-1);
    num_activated_T_pertime(n) = num_activated_T_pertime(n-1);
    T_dis2center_track_3severalT(:,n) = T_dis2center_track_3severalT(:,n-1);
    if index == 2 
        if num_tumor_cell_real <=0
            break
        else
            max_divcell_selection_attempts = 80;
            max_div_position_attempts = 200;
            num_newcell_per_div = 1;
            tumor_counts_per_clone(:,n) = tumor_counts_per_clone(:,n-1);
            num_tumor_cell_pertime(n) = num_tumor_cell_pertime(n-1);
            T_counts_per_clone(:,n) = T_counts_per_clone(:,n-1);
            T_dis2center_track(n) = T_dis2center_track(n-1);
            success_div = false;
            no_dividablecell = false;
            for div_cell_selection_attempt = 1:max_divcell_selection_attempts
                if ~success_div
                    div_clone_prob = nan(1,length(clone_group_ID));
                    num_cells_per_clone = cellfun(@length,clone_group_ID);
                    for i = 1:length(clone_group_ID)
                        if ~isempty(clone_group_ID{i})
                            div_clone_prob(i) = (averaged_tumor_div_rate_perclone(i)*length(clone_group_ID{i}))./sum((averaged_tumor_div_rate_perclone.*num_cells_per_clone));
                        else
                            div_clone_prob(i) = 0;
                        end
                    end
                    index_div_clone = [];
                    valid_tumor_idx = [];
                    tried_div_clone_idx = [];
                    while isempty(valid_tumor_idx)
                        if length(tried_div_clone_idx) == length(clone_group_ID)
                            no_dividablecell = true;
                            break
                        end
                        
                        nonempty_clone_id = find(~cellfun(@(x) all(isempty(x)), clone_group_ID));
                        index_div_clone = randsample(nonempty_clone_id, 1, true, div_clone_prob(nonempty_clone_id));
                        if ismember(index_div_clone, tried_div_clone_idx)
                            continue; 
                        end
                        tried_div_clone_idx(end + 1) = index_div_clone;
                        valid_tumor_idx = setdiff(intersect(clone_group_ID{index_div_clone},dividable_tumor_idx),killed_tumor_idx);
                    end
                    
                    
                    if ~no_dividablecell
                        index_div_cell = valid_tumor_idx(randsample(1:length(valid_tumor_idx),1));
                        
                        for div_attempt = 1: max_div_position_attempts
                            if ~success_div
                                div_angle = 2 * pi * rand(1,num_newcell_per_div);
                                new_divcell_x = tumor_positionhis_x(index_div_cell) + 2 * r_tumor * cos(div_angle);
                                new_divcell_y = tumor_positionhis_y(index_div_cell) + 2 * r_tumor * sin(div_angle);
                                new_divcell_position = [new_divcell_x, new_divcell_y];
                                newdiv_2center_dis = norm(new_divcell_position - center);
                                
                                if newdiv_2center_dis >= r_tumor_core + r_tumor
                                    newdiv2existingcell_dis = sqrt((new_divcell_x - tumor_positionhis_x).^2 + (new_divcell_y - tumor_positionhis_y).^2);
                                    if all(newdiv2existingcell_dis(~isnan(newdiv2existingcell_dis)) >= 2 * r_tumor)
                                         success_div = true;
                                        break
                                    end
                                end
                            end
                            if div_attempt == max_div_position_attempts && success_div == false
                                non_dividable_tumor_idx = [non_dividable_tumor_idx, [index_div_cell,killed_tumor_idx]];
                                non_dividable_tumor_idx = unique(non_dividable_tumor_idx);
                            end
                        end
                    else
                        break
                    end
                else
                    break
                end
            end
            if isempty(killed_tumor_idx)
                dividable_tumor_idx = setdiff((1:num_tumor_cell_track),non_dividable_tumor_idx);
            else
                dividable_tumor_idx = setdiff(1:num_tumor_cell_track, [non_dividable_tumor_idx, killed_tumor_idx]);
            end
            
            if success_div
                
                num_tumor_cell_track = num_tumor_cell_track + num_newcell_per_div;
                num_tumor_cell_real = num_tumor_cell_real +1;
                num_tumor_cell_pertime(n) = num_tumor_cell_pertime(n) + 1;
                dividable_tumor_idx = [dividable_tumor_idx,num_tumor_cell_track];
                tumor_positionhis_x (num_tumor_cell_track) = new_divcell_x;
                tumor_positionhis_y (num_tumor_cell_track) = new_divcell_y;
                tumor_div_rate_percell(end+1) = tumor_div_rate;
                tumor_mov_rate_percell(end+1) = tumor_mov_rate;     
                parent_cell_clone = find(cellfun(@(x) any(x == index_div_cell), clone_group_ID));
                parent_cell_antigen = clonegroup_antigen(parent_cell_clone);
                parent_cell_antigen = cell2mat(parent_cell_antigen);
                daughter_parent(end+1,:) = [num_tumor_cell_track,index_div_cell];
            end
        end
    end
end