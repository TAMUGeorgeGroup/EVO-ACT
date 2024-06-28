%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% author: Vicky Fan
% The following code is for tumor cell migration.
% All the values shown below are arbitrary examples. Please feel free to change them according to needs.
% Texas A&M University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if num_tumor_cell_real >0
    max_movetumor_selection_attempts = 50;
    tumor_move_dis = 2 * D * sqrt(tumor_mov_waiting);
    tumor_search_radius = 1.5;
    tumor_counts_per_clone(:,n) = tumor_counts_per_clone(:,n-1);
    num_tumor_cell_pertime(n) = num_tumor_cell_pertime(n-1);
    T_counts_per_clone(:,n) = T_counts_per_clone(:,n-1);
    T_dis2center_track(n) = T_dis2center_track(n-1);
    
    success_tumor_move = false(max_movetumor_selection_attempts,1);
    index_move_tumor_values = zeros(max_movetumor_selection_attempts,1);
    new_move_tumorcell_x_values = zeros(max_movetumor_selection_attempts,1);
    new_move_tumorcell_y_values = zeros(max_movetumor_selection_attempts,1);
    fiber_midpoint_x_local = fiber_midpoint_matrix(:,1);
    fiber_midpoint_y_local = fiber_midpoint_matrix(:,2);
    
    for move_tumor_selection_attempt = 1:max_movetumor_selection_attempts
        mean_direction=[];
        if ~success_tumor_move(move_tumor_selection_attempt)
           
            live_tumor_idx = find(~isnan(tumor_positionhis_x));
            index_move_tumor = live_tumor_idx(randi(length(live_tumor_idx)));
            while isnan(tumor_positionhis_x(index_move_tumor)) 
                index_move_tumor = randi(num_tumor_cell_track);
            end
            
            for tumor_move_attempt = 1:max_tumor_moveattempts
                dis_fiber2tumor = sqrt((fiber_midpoint_x_local-tumor_positionhis_x(index_move_tumor)).^2 + (fiber_midpoint_y_local-tumor_positionhis_y(index_move_tumor)).^2);
                fiber_within_radius_idx = find(dis_fiber2tumor <= tumor_search_radius);
                if isempty(fiber_within_radius_idx)
                    mean_direction = 2*pi*rand;
                else
                    mean_direction = fiber_direction_matrix(fiber_within_radius_idx(randi(length(fiber_within_radius_idx))));
                    if rand <0.5 
                        mean_direction = mod(mean_direction + pi,2*pi);
                    end
                end
                
                
                sigma_move = 0.01;
                untruncated_move_dir = makedist('Normal', mean_direction, sigma_move);
                truncated_mov_dir = truncate(untruncated_move_dir,0,2*pi);
                move_angle = random(truncated_mov_dir, 1);
                new_move_tumorcell_x = tumor_positionhis_x(index_move_tumor) + tumor_move_dis * cos(move_angle);
                new_move_tumorcell_y = tumor_positionhis_y(index_move_tumor) + tumor_move_dis * sin(move_angle);
                
                num_overlapcheck_en_route = 5;
                
                dis_tumor2tumor = sqrt((new_move_tumorcell_x - tumor_positionhis_x).^2 + (new_move_tumorcell_y - tumor_positionhis_y).^2);
                condition = all(dis_tumor2tumor > tumor_move_dis, 1);
                condition(index_move_tumor) = false;
                condition_idx = find(condition);
                row_indices = find(all(dis_tumor2tumor(condition_idx) > tumor_move_dis, 1));
                num_tumor_within_movingradius = length(condition_idx);
                
                if num_tumor_within_movingradius >0
                    tumor_x_within_moveradius = tumor_positionhis_x(row_indices(1:num_tumor_within_movingradius));
                    tumor_y_within_moveradius = tumor_positionhis_y(row_indices(1:num_tumor_within_movingradius));
                else
                    tumor_x_within_moveradius = [];
                    tumor_y_within_moveradius = [];
                end
                
                newmovetumor2center_dis = nan(num_overlapcheck_en_route,1);
                newmovetumor2existingcell_dis = {};
                
                for check_step = 1:num_overlapcheck_en_route
                    m = check_step / num_overlapcheck_en_route;
                    intermediate_x_tumor = (1 - m) * tumor_positionhis_x(index_move_tumor) + m * new_move_tumorcell_x;
                    intermediate_y_tumor = (1 - m) * tumor_positionhis_y(index_move_tumor) + m * new_move_tumorcell_y;
                    newmovetumor2center_dis(check_step) = sqrt((intermediate_x_tumor).^2 + (intermediate_y_tumor).^2);
                    
                    if num_tumor_within_movingradius == 0
                        condition1 = true; 
                    else
                        newmovetumor2existingcell_dis{check_step} = sqrt((intermediate_x_tumor - tumor_x_within_moveradius).^2 + (intermediate_y_tumor - tumor_y_within_moveradius).^2); 
                        newmovetumor2existingcell_dis_matrix = cell2mat(newmovetumor2existingcell_dis);
                        condition1 = all(newmovetumor2existingcell_dis_matrix > 2*r_tumor);
                    end
                end
                condition2 = all(newmovetumor2center_dis >= r_tumor_core + r_tumor);
                
                if condition1 && condition2
                    success_tumor_move(move_tumor_selection_attempt) = true;
                    index_move_tumor_values(move_tumor_selection_attempt) = index_move_tumor;
                    new_move_tumorcell_x_values(move_tumor_selection_attempt) = new_move_tumorcell_x;
                    new_move_tumorcell_y_values(move_tumor_selection_attempt) = new_move_tumorcell_y;
                    break;
                end
            end
        end
    end
    
    success_tumor_move_attempt = find(success_tumor_move ~= 0, 1, 'first');
    
    if success_tumor_move(success_tumor_move_attempt) 
        final_move_tumor_idx = index_move_tumor_values(success_tumor_move_attempt);
        final_tumor_x = new_move_tumorcell_x_values(success_tumor_move_attempt);
        final_tumor_y = new_move_tumorcell_y_values(success_tumor_move_attempt);
        tumor_positionhis_x(final_move_tumor_idx) = final_tumor_x;
        tumor_positionhis_y(final_move_tumor_idx) = final_tumor_y;
    end
    
end
