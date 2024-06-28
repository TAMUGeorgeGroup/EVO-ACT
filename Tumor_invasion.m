%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% author: Vicky Fan
% The following code is for tumor invasion.
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
    if index == 1 % Tumor cell invades.
        if num_tumor_cell_real <= num_tumor_pheri
            parfor inv_attempt = 1: max_invattempts
                if ~success_inv(inv_attempt)
                    inv_angle = 2 * pi * rand;
                    new_invcell_x = (r_tumor_core + r_tumor) * cos(inv_angle);
                    new_invcell_y = (r_tumor_core + r_tumor) * sin(inv_angle);
                    newinv2existingcell_dis = sqrt((new_invcell_x - tumor_positionhis_x_local).^2 + (new_invcell_y - tumor_positionhis_y_local).^2);
                    if all(newinv2existingcell_dis >= 2*r_tumor)
                        success_inv(inv_attempt) = true;
                        new_invcell_x_values(inv_attempt) = new_invcell_x;
                        new_invcell_y_values(inv_attempt) = new_invcell_y;
                    end
                end
            end
        end
        success_inv_attempt = find(success_inv ~= 0, 1, 'first');
        if success_inv(success_inv_attempt) 
            final_new_invcell_x = new_invcell_x_values(success_inv_attempt);
            final_new_invcell_y = new_invcell_y_values(success_inv_attempt);
            num_tumor_cell_track = num_tumor_cell_track + 1;
            num_tumor_cell_real = num_tumor_cell_real + 1;
            tumor_positionhis_x (num_tumor_cell_track, :) = nan;
            tumor_positionhis_x (num_tumor_cell_track, n) = final_new_invcell_x;
            tumor_positionhis_y (num_tumor_cell_track, :) = nan;
            tumor_positionhis_y (num_tumor_cell_track, n) = final_new_invcell_y;
            clone_group_ID{1} = [clone_group_ID{1}, num_tumor_cell_track];
        end
    end
    
end