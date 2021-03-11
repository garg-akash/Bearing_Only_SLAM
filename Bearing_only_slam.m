close all
clear
clc

#load dependencies
addpath "./g2o_wrapper"
source './common_stuff.m'
source './bearing_least_squares_pose.m'
source './bearing_least_squares_landmark.m'
source './bearing_least_squares_solve.m'

#load ground truth dataset
[landmarks, poses, transitions, observations] = loadG2o('slam2D_bearing_only_ground_truth.g2o');

num_landmarks = size(landmarks,2);
num_poses = size(poses,2);
num_transitions = size(transitions,2);
num_observations = size(observations,2);

XR_true = zeros(3,3,num_poses);
XL_true = zeros(2,num_landmarks);

state_to_id_xr = ones(num_poses, 1)*-1;
state_to_id_xl = ones(num_landmarks, 1)*-1;

global pose_dim = 3;
global landmark_dim = 2;

for i=1:num_poses
  XR_true(:,:,i) = v2t([poses(i).x,poses(i).y,poses(i).theta]);
  state_to_id_xr(i) = poses(i).id;
endfor

for i=1:num_landmarks
  XL_true(:,i) = [landmarks(i).x_pose;landmarks(i).y_pose];
  state_to_id_xl(i) = landmarks(i).id;
endfor

#load initial guess dataset
[landmarks, poses, transitions, observations] = loadG2o('slam2D_bearing_only_initial_guess.g2o');
XR_guess = zeros(3,3,num_poses);
for i=1:num_poses
	XR_guess(:,:,i)=v2t([poses(i).x poses(i).y poses(i).theta]);
	state_to_id_xr(i) = poses(i).id;
endfor;
#count total landmark observations
num_landmark_measurements = 0;
for i=1:num_observations
  num_landmark_measurements += size(observations(i).observation,2);
endfor
#robot-landmark observation
Zl = zeros(1,num_landmark_measurements);
landmark_associations = zeros(2,num_landmark_measurements); #(observing robot position,observed landmark position) 

measurement_num=1;
for i=1:num_observations
  pose_state = id_to_state(state_to_id_xr,observations(i).pose_id);
  for j=1:size(observations(i).observation,2)
    landmark_state = id_to_state(state_to_id_xl,observations(i).observation(j).id);
	  Zl(1,measurement_num) = observations(i).observation(j).bearing;
    landmark_associations(:,measurement_num) = [pose_state;landmark_state];
    measurement_num++;
  endfor;
endfor;

#count total transitions
num_transitions = size(transitions,2);

#robot-robot observations
Zr = zeros(3,3,num_transitions);
pose_associations = zeros(2,num_transitions);

for i=1:num_transitions
  Zr(:,:,i) = v2t(transitions(i).v);
  pose_state_from = id_to_state(state_to_id_xr,transitions(i).id_from);
  pose_state_to = id_to_state(state_to_id_xr,transitions(i).id_to);
  pose_associations(:,i) = [pose_state_from;pose_state_to];
endfor;

XL_guess = zeros(2,num_landmarks);
#Triangulation starts
[sort_landmark_state sort_landmark_index] = sort(landmark_associations(2,:));
landmark_unique = unique(sort_landmark_state);
landmark_unique_count = histc(sort_landmark_state,landmark_unique);
for i=1:length(landmark_unique_count)
  r_vec = zeros(3,landmark_unique_count(i));
  landmark_unique_measure = zeros(1,landmark_unique_count(i));
  for j=1:landmark_unique_count(i)
    r_state = landmark_associations(1,sort_landmark_index(sum(landmark_unique_count(1:i))-landmark_unique_count(i)+j));
    r_vec(:,j) = t2v(XR_guess(:,:,r_state));
    landmark_unique_measure(1,j) = Zl(1,sort_landmark_index(sum(landmark_unique_count(1:i))-landmark_unique_count(i)+j));
    r_vec(3,j) = r_vec(3,j) + landmark_unique_measure(1,j);
    r_vec(3,j) = atan2(sin(r_vec(3,j)),cos(r_vec(3,j)));
  endfor;
  if (landmark_unique_count(i)==1)  
    XL_guess(:,i) = [r_vec(1,1)+cos(r_vec(3,1));r_vec(2,1)+sin(r_vec(3,1))];  
  else 
    heading_diff = r_vec(3,:)-r_vec(3,:)';
    heading_diff = atan2(sin(heading_diff),cos(heading_diff));    
    heading_diff = sin(abs(heading_diff));
    [rows,cols] = find(heading_diff==max(heading_diff(:))); #find two poses having max difference in their measurement heading
    row = rows(1); col = cols(1);
    start_pt1 = r_vec(1:2,row);
    end_pt1 = r_vec(1:2,row) + [cos(r_vec(3,row));sin(r_vec(3,row))];
    start_pt2 = r_vec(1:2,col);
    end_pt2 = r_vec(1:2,col) + [cos(r_vec(3,col));sin(r_vec(3,col))];
    XL_guess(:,i) = computeIntersection(start_pt1,end_pt1,start_pt2,end_pt2);
  endif;
endfor;
#Triangulation ends

%XR_guess = XR_true;
%XL_guess = XL_true;

#Call Solver
damping=1e-4;
kernel_threshold=10;
num_iterations=20;
[XR, XL, chi_stats_l, num_inliers_l, chi_stats_r, num_inliers_r, H, b] = doTotalLS(XR_guess,XL_guess,Zl,landmark_associations,
                                                                            Zr,pose_associations,num_iterations,damping,kernel_threshold);

# Plot State
figure(1);
hold on;
grid;

subplot(2,2,1);
title("Landmark Initial Guess");
plot(XL_true(1,:),XL_true(2,:),'b*',"linewidth",2);
hold on;
plot(XL_guess(1,:),XL_guess(2,:),'ro',"linewidth",2);
legend("Landmark True", "Guess");grid;


subplot(2,2,2);
title("Landmark After Optimization");
plot(XL_true(1,:),XL_true(2,:),'b*',"linewidth",2);
hold on;
plot(XL(1,:),XL(2,:),'ro',"linewidth",2);
legend("Landmark True", "Guess");grid;


subplot(2,2,3);
title("Poses Initial Guess");
plot(XR_true(1,:),XR_true(2,:),'b*',"linewidth",2);
hold on;
plot(XR_guess(1,:),XR_guess(2,:),'ro',"linewidth",2);
legend("Poses True", "Guess");grid;


subplot(2,2,4);
title("Poses After Optimization");
plot(XR_true(1,:),XR_true(2,:),'b*',"linewidth",2);
hold on;
plot(XR(1,:),XR(2,:),'ro',"linewidth",2);
legend("Poses True", "Guess"); grid;


figure(2);
hold on;
grid;
title("chi evolution");

subplot(2,2,1);
plot(chi_stats_r, 'r-', "linewidth", 2);
legend("Chi Poses"); grid; xlabel("iterations");
subplot(2,2,2);
plot(num_inliers_r, 'b-', "linewidth", 2);
legend("#inliers"); grid; xlabel("iterations");

subplot(2,2,3);
plot(chi_stats_l, 'r-', "linewidth", 2);
legend("Chi Landmark"); grid; xlabel("iterations");
subplot(2,2,4);
plot(num_inliers_l, 'b-', "linewidth", 2);
legend("#inliers"); grid; xlabel("iterations");

figure(3);
title("H matrix");
H_ =  H./H;                      # NaN and 1 element
H_(isnan(H_))=0;                 # Nan to Zero
H_ = abs(ones(size(H_)) - H_);   # switch zero and one
H_ = flipud(H_);                 # switch rows
colormap(gray(64));
hold on;
image([0.5, size(H_,2)-0.5], [0.5, size(H_,1)-0.5], H_*64);
hold off;