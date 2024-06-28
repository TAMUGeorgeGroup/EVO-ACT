%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% author: Vicky Fan
% The following code is for T cell migration.
% All the values shown below are arbitury examples. Please feel free to change them by needs.
% Texas A&M University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for move_T_selection_attempt = 1:max_moveT_selection_attempts
    if move_T_selection_attempt >1
        if success_T_move(move_T_selection_attempt-1)
            break
        end
    end
    
    index_move_T = live_T_idx(randi(length(live_T_idx)));
    
    for T_move_attempt = 1:max_T_moveattempts
        dis2T = sqrt((fiber_midpoint_matrix(:,1)-T_positionhis_x(index_move_T)).^2 + (fiber_midpoint_matrix(:,2)-T_positionhis_y(index_move_T)).^2);
        fiber_check_direction_T = find(dis2T <= T_search_radius);
        if isempty(fiber_check_direction_T)
            mean_direction_T = 2*pi*rand;
        else
            fiber_chosen_idx = fiber_check_direction_T(randi(length(fiber_check_direction_T)));
            mean_direction_T = fiber_direction_matrix(fiber_chosen_idx);
            if rand <0.5
                mean_direction_T = mod(mean_direction_T + pi,2*pi);
            end
        end
        
        gradient_sum = [0, 0];
        fiber_search_radius = T_search_radius;
        num_fibers_within_radius = length(fiber_check_direction_T);
        D_T = 50000;
        nonNan_tumor_idx = find(~isnan(tumor_positionhis_x));
        
        if all(isnan(tumor_positionhis_x))
            T_move_dis = D_T * sqrt(T_move_waiting);
            T_move_direction = 2*pi*rand;
            new_move_Tcell_x = T_positionhis_x(index_move_T) + T_move_dis * cos(T_move_direction);
            new_move_Tcell_y = T_positionhis_y(index_move_T) + T_move_dis * sin(T_move_direction);
        else
            d = sqrt((T_positionhis_x(index_move_T) - tumor_positionhis_x(nonNan_tumor_idx)).^2 + (T_positionhis_y(index_move_T) - tumor_positionhis_y(nonNan_tumor_idx)).^2);
            gradient_magnitude = -80*exp(-d/110);
            vector = [T_positionhis_x(index_move_T) - tumor_positionhis_x(nonNan_tumor_idx) T_positionhis_y(index_move_T) - tumor_positionhis_y(nonNan_tumor_idx)];
            gradient_sum = gradient_magnitude .* vector / norm(vector);
            if num_fibers_within_radius > 0
                T_move_dis = 300*norm(gradient_sum) * T_move_waiting/sqrt(num_fibers_within_radius);
            else
                T_move_dis = 300*norm(gradient_sum) * T_move_waiting;
            end
            average_tumor_position = [sum(tumor_positionhis_x(nonNan_tumor_idx))/length(nonNan_tumor_idx), sum(tumor_positionhis_y(nonNan_tumor_idx))/length(nonNan_tumor_idx)];
            T_move_direction_1 = mean_direction_T;
            new_move_Tcell_x_1 = T_positionhis_x(index_move_T) + T_move_dis * cos(T_move_direction_1);
            new_move_Tcell_y_1 = T_positionhis_y(index_move_T) + T_move_dis * sin(T_move_direction_1);
            new_move_T_position_1 = [new_move_Tcell_x_1, new_move_Tcell_y_1];
            
            T_move_direction_2 = mod(mean_direction_T - pi, 2*pi);
            new_move_Tcell_x_2 = T_positionhis_x(index_move_T) + T_move_dis * cos(T_move_direction_2);
            new_move_Tcell_y_2 = T_positionhis_y(index_move_T) + T_move_dis * sin(T_move_direction_2);
            new_move_T_position_2 = [new_move_Tcell_x_2, new_move_Tcell_y_2];
            
            if norm(new_move_T_position_1 - average_tumor_position) < norm(new_move_T_position_2 - average_tumor_position)
                T_move_direction = T_move_direction_1;
                new_move_Tcell_x = new_move_Tcell_x_1;
                new_move_Tcell_y = new_move_Tcell_y_1;
            elseif norm(new_move_T_position_1 - average_tumor_position) >= norm(new_move_T_position_2 - average_tumor_position)
                T_move_direction = T_move_direction_2;
                new_move_Tcell_x = new_move_Tcell_x_2;
                new_move_Tcell_y = new_move_Tcell_y_2;
            end
        end
        for check_step = 1:num_ToverlapT_check_en_route
            m_T = check_step / num_ToverlapT_check_en_route;
            intermediate_x_T = (1- m_T) * T_positionhis_x(index_move_T) + m_T * new_move_Tcell_x;
            intermediate_y_T = (1- m_T) * T_positionhis_y(index_move_T) + m_T * new_move_Tcell_y;
            newmoveT_2existingT_dis = sqrt((intermediate_x_T - T_positionhis_x).^2 + (intermediate_y_T - T_positionhis_y).^2);
            if all(newmoveT_2existingT_dis) >= 2*r_T
                T_move(check_step) = true;
            end
        end
        
        if all(T_move) == true
            success_T_move(move_T_selection_attempt) = true;
            index_move_T_values(move_T_selection_attempt) = index_move_T;
            break
        end
        
    end
end