%% TACS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% author: Vicky Fan
% The following code is used to generate a randomly packed TACS1 fiber network. 
% Texas A&M University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


center = [0,0];
r_matrix = 40;
num_fibers_inner = 5000;
mu_length=0.2;
sigma_length=0.1;
r_tumor_core = 5;
mu_beta = 0;
R_max = r_matrix;

b = 2*num_fibers_inner/(R_max - r_tumor_core);
a = b/(R_max - r_tumor_core);

fiber_midpoint_matrix = [];
fiber_beta_matrix = [];
fiber_start_matrix = [];
fiber_end_matrix = [];
fiber_direction_matrix = [];
fiber_idx = 1;
increment = 0.1;
num_fiber_per_circle = nan(1,(R_max-r_tumor_core)/increment);

% figure;
% rectangle('Position',[center(1) - r_matrix, center(2) - r_matrix, 2*r_matrix, 2*r_matrix],'Curvature',[1,1],'FaceColor',[0.8, 0.8, 0.8]);
% viscircles(center,r_matrix,'color','w','LineWidth',1);
% set(gca, 'YDir', 'normal');
% hold on

for R=r_tumor_core:increment:R_max
    circle_id = round((R - r_tumor_core + increment)/increment);
    num_fibers_atR = round(-a*(R-r_tumor_core)+b);
    num_fiber_per_circle(circle_id) = num_fibers_atR;
    fiber_length = normrnd(mu_length,sigma_length,[1,num_fibers_atR]);
    alpha_list = 2*pi*rand(1,num_fibers_atR);
    for i = 1:num_fibers_atR
        alpha = alpha_list(i);
        fiber_midpoint = [R*cos(alpha) R*sin(alpha)];
        fiber_direction = 2*pi*rand; 
        fiber_start = fiber_midpoint + (fiber_length(i)/2) .* [cos(fiber_direction), sin(fiber_direction)];
        fiber_end = fiber_midpoint - (fiber_length(i)/2) .* [cos(fiber_direction), sin(fiber_direction)];
        fibermid_2center_dis = norm(fiber_midpoint);
        fiberstart_2center_dis = norm(fiber_start);
        fiberend_2center_dis = norm(fiber_end);
        fiber2_center_dis_vector = [fibermid_2center_dis,fiberstart_2center_dis,fiberend_2center_dis];
        fiber_midpoint_matrix(fiber_idx,:) = fiber_midpoint;
        fiber_start_matrix(fiber_idx,:) = fiber_start;
        fiber_end_matrix(fiber_idx,:) = fiber_end;
        fiber_direction_matrix(fiber_idx) = fiber_direction;
        fiber_length_matrix(fiber_idx) = fiber_length(i);
        fiber_idx = fiber_idx +1;
    end
end

%% 
%%