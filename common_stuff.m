1;
function R=rotation2D(theta)
  s = sin(theta);
  c = cos(theta);
  R = [c -s;
       s  c];
endfunction;

function T=v2t(v)
  T = eye(3);
  T(1:2,1:2) = rotation2D(v(3));
  T(1:2,3) = v(1:2);
endfunction;

function v=t2v(A)
  v(1:2, 1) = A(1:2,3);
  v(3,1) = atan2(A(2,1),A(1,1));
endfunction;

function state=id_to_state(query_vector,query_id)
  state = find(query_vector==query_id);
  if (size(state,2)==0)
    state = -1;
  endif
endfunction;

# retrieves the index in the perturbation vector, that corresponds to
# a certain pose
# input:
#   pose_index:     the index of the pose for which we want to compute the
#                   index
#   num_poses:      number of pose variables in the state
#   num_landmarks:  number of pose variables in the state
# output:
#   v_idx: the index of the sub-vector corrsponding to 
#          pose_index, in the array of perturbations  (-1 if error)
function v_idx=poseMatrixIndex(pose_index, num_poses, num_landmarks)
  global pose_dim;
  global landmark_dim;
  if (pose_index>num_poses)
    v_idx = -1;
    return;
  endif;
  v_idx = 1+(pose_index-1)*pose_dim;
endfunction;

# retrieves the index in the perturbation vector, that corresponds to
# a certain landmark
# input:
#   landmark_index:     the index of the landmark for which we want to compute the
#                   index
#   num_poses:      number of pose variables in the state
#   num_landmarks:  number of pose variables in the state
# output:
#   v_idx: the index of the perturnation corrsponding to the
#           landmark_index, in the array of perturbations

function v_idx=landmarkMatrixIndex(landmark_index, num_poses, num_landmarks)
  global pose_dim;
  global landmark_dim;
  if (landmark_index>num_landmarks)
    v_idx = -1;
    return;
  endif;
  v_idx = 1 + (num_poses)*pose_dim + (landmark_index-1) * landmark_dim;
endfunction;

function v=flattenIsometryByColumns(T)
  v=zeros(6,1);
  v(1:4)=reshape(T(1:2,1:2),4,1);
  v(5:6)=T(1:2,3);
endfunction;

function T=unflattenIsometryByColumns(v)
  T=eye(3);
  T(1:2,1:2)=reshape(v(1:4),2,2);
  T(1:2,3)=v(5:6);
endfunction;

#pt1 and pt2 are the starting and ending points of line 1
#pt3 and pt4 are the starting and ending points of line 2
function intersection_pt = computeIntersection(pt1,pt2,pt3,pt4)
  m1 = (pt2(2)-pt1(2))/(pt2(1)-pt1(1));
  c1 = (pt2(1)*pt1(2)-pt1(1)*pt2(2))/(pt2(1)-pt1(1));
  
  m2 = (pt4(2)-pt3(2))/(pt4(1)-pt3(1));
  c2 = (pt4(1)*pt3(2)-pt3(1)*pt4(2))/(pt4(1)-pt3(1));
  
  x = (c2-c1)/(m1-m2);
  y = m1*x + c1;
  intersection_pt = [x;y];
endfunction;